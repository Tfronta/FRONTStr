"""Orchestrator: evidence clusters → :class:`MarkerResult`.

This is the public entry point of the Interpretation layer. Higher-level
code (CLI, report, exports) calls :func:`interpret_marker` (single locus)
or :func:`interpret_run` (full panel).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
from pathlib import Path

from frontstr.evidence.cluster import Cluster, cluster_observations
from frontstr.evidence.consensus import poa_backend_name
from frontstr.evidence.pileup import _DEFAULT_FLANK_ANCHOR, PileupCounts, pileup_locus
from frontstr.interp.allele_numeric import compute_allele_numeric, resolve_ref_anchor_bp
from frontstr.interp.amel import interpret_amel
from frontstr.interp.catalog import annotate_alleles
from frontstr.interp.classify import classify_allele
from frontstr.interp.flags import derive_marker_flags
from frontstr.interp.haplotype import suppress_hp_phantoms
from frontstr.interp.isfg import ce_from_brackets, ce_from_length, compress_isfg
from frontstr.interp.models import Allele, CallRule, Flag, FlagCode, MarkerResult, TriType
from frontstr.interp.naming import NameStatus, StrNamer, default_namer
from frontstr.interp.qc import QcThresholds, derive_run_qc_flags
from frontstr.interp.stutter import build_expected_stutter
from frontstr.interp.triallelic import DEFAULT_MIN_PHR_FOR_HET, call_profile
from frontstr.log import get_logger
from frontstr.motifs import repeat_core_span, reverse_complement
from frontstr.panel.catalog import AlleleCatalog
from frontstr.panel.models import Panel, System
from frontstr.params import RunParameters
from frontstr.trace import FLANK_SHOWN, BinTrace, ClusterTrace, LocusTrace

log = get_logger(__name__)

DEFAULT_ANALYTICAL_THRESH = 0.02
DEFAULT_CALLING_THRESH = 0.10
DEFAULT_PARENT_FRACTION = 0.20


def interpret_marker(
    *,
    system: System,
    clusters: list[Cluster],
    analytical_thresh: float = DEFAULT_ANALYTICAL_THRESH,
    calling_thresh: float = DEFAULT_CALLING_THRESH,
    parent_fraction: float = DEFAULT_PARENT_FRACTION,
    ref_length_bp: int | None = None,
    catalog: AlleleCatalog | None = None,
    namer: StrNamer | None = None,
    trace: LocusTrace | None = None,
    min_phr_for_het: float = DEFAULT_MIN_PHR_FOR_HET,
    min_reads_third: int | None = None,
) -> MarkerResult:
    """Interpret one marker's evidence into a :class:`MarkerResult`.

    Args:
        system: Marker definition.
        clusters: Output of :func:`frontstr.evidence.cluster.cluster_observations`.
        analytical_thresh: Below this fraction → noise.
        calling_thresh: Below this fraction (but >= analytical) → artefact.
        parent_fraction: Clusters with fraction ≥ this become parents for
            stutter-expectation calculation (default 0.20).
        ref_length_bp: Optional override for reference TR length (bp) used for
            ``bp_diff`` and compound numeric alleles. If ``None``, uses the
            panel span (``ref_end - ref_start + 1``).
        namer: STRNaming namer supplying the canonical allele number. Defaults
            to the bundled one. Markers it has no reporting range for, and
            consensuses it cannot locate the range in, fall back to the legacy
            CE per allele — so this never needs disabling for coverage reasons.
        trace: Optional :class:`frontstr.trace.LocusTrace` to fill with the
            interpretation half of the per-locus narrative. The evidence half is
            filled by :func:`_safe_pileup_and_cluster`.
        min_phr_for_het: Minor allele as a fraction of the major before a
            heterozygote is called on read counts alone.
        min_reads_third: Absolute read floor for a 3rd candidate; ``None`` uses
            each marker's own ``min_reads_third``.
    """
    if namer is None:
        namer = default_namer()
    total_reads = sum(c.n_reads for c in clusters)
    ref_length_bp = resolve_ref_anchor_bp(system, explicit=ref_length_bp)

    alleles = [
        _allele_from_cluster(idx, c, system, ref_length_bp, namer) for idx, c in enumerate(clusters)
    ]

    parents = [
        c
        for c, a in zip(clusters, alleles, strict=True)
        if a.fraction(total_reads) >= parent_fraction and not a.is_deletion
    ]
    expected = build_expected_stutter(parents, system)

    for a in alleles:
        a.expected_stutter = expected.get(a.consensus, 0.0)
        a.status = classify_allele(
            a,
            total_reads=total_reads,
            expected_stutter=expected,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
        )

    # Suppress same-haplotype split-allele phantoms before the profile is
    # called, so a phantom can neither be reported nor raise a false mixture.
    # No-op on unphased BAMs and on allow_triallelic markers.
    n_phantoms = suppress_hp_phantoms(alleles, system)

    # Enrich ISFG / CE / iso-allele suffix from the curated catalog (no-op if None).
    annotate_alleles(alleles, system, catalog)

    alleles_called, call_rule, tri_type = call_profile(
        alleles,
        system,
        calling_thresh=calling_thresh,
        min_phr_for_het=min_phr_for_het,
        min_reads_third=min_reads_third,
    )
    _report_naming_fallbacks(system, alleles, alleles_called)
    if trace is not None:
        _fill_interp_trace(
            trace, system, alleles, alleles_called, total_reads, n_phantoms, call_rule, tri_type
        )

    result = MarkerResult(
        marker_name=system.name,
        system=system,
        alleles=alleles,
        alleles_called=alleles_called,
        call_rule=call_rule,
        tri_type=tri_type,
        total_reads=total_reads,
        expected_stutter=expected,
        analytical_thresh=analytical_thresh,
        calling_thresh=calling_thresh,
    )
    derive_marker_flags(result)
    if trace is not None:
        trace.flags = [(f.severity.value, f.code.value) for f in result.flags]
        trace.allele_balance = result.allele_balance
    return result


def interpret_run(
    *,
    bam: Path,
    panel: Panel,
    min_mapq: int = 20,
    identity_threshold: float = 0.97,
    len_tolerance_bp: int = 0,
    analytical_thresh: float = DEFAULT_ANALYTICAL_THRESH,
    calling_thresh: float = DEFAULT_CALLING_THRESH,
    reference_fasta: Path | None = None,
    catalog: AlleleCatalog | None = None,
    qc_thresholds: QcThresholds | None = None,
    on_trace: Callable[[LocusTrace], None] | None = None,
    params: RunParameters | None = None,
) -> list[MarkerResult]:
    """End-to-end: for each marker in ``panel``, pileup → cluster → interpret.

    Args:
        bam: Indexed sample BAM or CRAM.
        panel: Panel definition.
        reference_fasta: Reference FASTA path; required when ``bam`` is a CRAM.
        qc_thresholds: Laboratory QC policy for the run-level flags (coverage,
            strand bias). Defaults apply when omitted.
        on_trace: Called with a filled :class:`frontstr.trace.LocusTrace` as each
            locus finishes, so a CLI can narrate the run as it happens rather
            than after it. ``None`` skips trace collection entirely.
        params: The run's :class:`frontstr.params.RunParameters`. When given it
            **overrides** the individual keyword arguments above, so there is one
            source of truth for what a run used; the keywords remain for library
            callers that only need one or two.

    Returns:
        One :class:`MarkerResult` per marker in panel order, each carrying its
        marker-level *and* run-level QC flags. Empty pileups produce a
        ``NO_DATA`` result rather than failing.
    """
    out: list[MarkerResult] = []
    if params is not None:
        min_mapq = params["min_mapq"]
        identity_threshold = params["identity_threshold"]
        len_tolerance_bp = params["len_tolerance_bp"]
        analytical_thresh = params["analytical_thresh"]
        calling_thresh = params["calling_thresh"]
    min_phr_for_het = params["min_phr_for_het"] if params else DEFAULT_MIN_PHR_FOR_HET
    min_reads_third = params["min_reads_third"] if params else None
    flank_anchor = params["flank_anchor"] if params else _DEFAULT_FLANK_ANCHOR
    # Built once per run: seeding the reference structures is the expensive part
    # and the result is immutable across markers.
    namer = default_namer()
    log.info(
        "run.start",
        bam=str(bam),
        panel=panel.name,
        panel_version=panel.version,
        n_markers=len(panel.systems),
        min_mapq=min_mapq,
        identity_threshold=identity_threshold,
        analytical_thresh=analytical_thresh,
        calling_thresh=calling_thresh,
        catalog=bool(catalog),
        strnaming=bool(namer),
    )
    for system in panel.systems:
        if system.marker_type == "amel":
            amel = interpret_amel(system, bam, min_mapq=min_mapq, reference_fasta=reference_fasta)
            out.append(amel)
            if on_trace is not None:
                # Traced too, so the run summary's locus count matches the panel.
                # Amelogenin is not a tandem repeat: no binning, no clustering,
                # no naming. Saying so is more useful than omitting it.
                amel_trace = LocusTrace(
                    marker=system.name,
                    chrom=system.chromosome,
                    start=system.ref_start,
                    end=system.ref_end,
                    motif=system.motif,
                    period=system.period,
                    strand=system.strand,
                    min_mapq=min_mapq,
                )
                amel_trace.note = (
                    "Sex typing, not a tandem repeat: counts reads at AMELX and "
                    "AMELY instead of binning and clustering (interp/amel.py)."
                )
                amel_trace.called_labels = [a.number_label for a in amel.alleles_called]
                amel_trace.call_rule = amel.call_rule.value
                amel_trace.flags = [(f.severity.value, f.code.value) for f in amel.flags]
                on_trace(amel_trace)
            continue
        trace = (
            LocusTrace(
                marker=system.name,
                chrom=system.chromosome,
                start=system.ref_start,
                end=system.ref_end,
                motif=system.motif,
                period=system.period,
                strand=system.strand,
                min_mapq=min_mapq,
                identity_threshold=identity_threshold,
                analytical_thresh=analytical_thresh,
                calling_thresh=calling_thresh,
            )
            if on_trace is not None
            else None
        )
        clusters = _safe_pileup_and_cluster(
            bam=bam,
            system=system,
            min_mapq=min_mapq,
            identity_threshold=identity_threshold,
            len_tolerance_bp=len_tolerance_bp,
            reference_fasta=reference_fasta,
            trace=trace,
            flank_anchor=flank_anchor,
        )
        result = interpret_marker(
            system=system,
            clusters=clusters,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
            catalog=catalog,
            namer=namer,
            trace=trace,
            min_phr_for_het=min_phr_for_het,
            min_reads_third=min_reads_third,
        )
        if on_trace is not None and trace is not None:
            on_trace(trace)
        log.debug(
            "marker.called",
            marker=system.name,
            call_rule=result.call_rule.value,
            total_reads=result.total_reads,
            n_clusters=len(clusters),
            alleles=[a.number_label for a in result.alleles_called],
            number_method=sorted({a.number_method for a in result.alleles_called}),
        )
        out.append(result)

    thresholds = derive_run_qc_flags(out, qc_thresholds)
    if params is not None:
        _flag_non_default_thresholds(out, params)
    log.info(
        "run.complete",
        n_markers=len(out),
        n_called=sum(1 for r in out if r.alleles_called),
        flags=dict(sorted(Counter(f.code.value for r in out for f in r.flags).items())),
        low_coverage_reads=thresholds.low_coverage_reads,
    )
    return out


def _split_for_display(consensus: str, system: System) -> tuple[str, str, str, bool]:
    """Split a consensus into ``(left flank, repeat core, right flank, found)``.

    Returned in canonical (motif) orientation so it reads the same way as the
    ISFG string — for a minus-strand marker the raw consensus is reference
    orientation and the motif does not appear in it at all.

    The core is whole and the flanks are trimmed to
    :data:`frontstr.trace.FLANK_SHOWN`: a panel window is roughly 80% flank, so
    printing all of it would bury the part that distinguishes two alleles.
    """
    if not consensus:
        return "", "", "", False
    seq = reverse_complement(consensus) if system.strand == "-" else consensus
    motifs = [m for m in system.motif.split(",") if m]
    span = repeat_core_span(seq, motifs) if motifs else None
    if span is None:
        return "", seq, "", False
    start, end = span
    left = seq[max(0, start - FLANK_SHOWN) : start]
    right = seq[end : end + FLANK_SHOWN]
    return (
        ("…" + left) if start > FLANK_SHOWN else left,
        seq[start:end],
        (right + "…") if end + FLANK_SHOWN < len(seq) else right,
        True,
    )


def _fill_interp_trace(
    trace: LocusTrace,
    system: System,
    alleles: list[Allele],
    called: list[Allele],
    total_reads: int,
    n_phantoms: int,
    call_rule: CallRule,
    tri_type: TriType,
) -> None:
    """Copy the interpretation layer's decisions onto the locus trace."""
    called_ids = {id(a) for a in called}
    splits = {a.cluster_index: _split_for_display(a.consensus, system) for a in alleles}
    trace.clusters = [
        ClusterTrace(
            index=a.cluster_index,
            n_reads=a.n_reads_total,
            fraction=a.fraction(total_reads),
            length_bp=a.length_bp,
            consensus_method=a.consensus_method,
            n_hp1=a.n_reads_hp1,
            n_hp2=a.n_reads_hp2,
            n_untagged=a.n_reads_hp_none,
            number_label=a.number_label,
            number_method=a.number_method,
            naming_status=a.strnaming_status,
            strnaming_name=a.strnaming_name,
            isfg=a.isfg,
            expected_stutter=a.expected_stutter,
            status=a.status.value,
            n_reads_absorbed=a.n_reads_absorbed,
            hp_rescued=a.hp_rescued,
            called=id(a) in called_ids,
            flank_left=splits[a.cluster_index][0],
            core=splits[a.cluster_index][1],
            flank_right=splits[a.cluster_index][2],
            core_found=splits[a.cluster_index][3],
            phase_set=a.phase_set,
            n_phase_sets=a.n_phase_sets,
        )
        for a in alleles
    ]
    trace.n_suppressed_phantoms = n_phantoms
    trace.called_labels = [a.number_label for a in called]
    trace.call_rule = call_rule.value
    trace.tri_type = tri_type.value


def _safe_pileup_and_cluster(
    *,
    bam: Path,
    system: System,
    min_mapq: int,
    identity_threshold: float,
    len_tolerance_bp: int,
    reference_fasta: Path | None = None,
    trace: LocusTrace | None = None,
    flank_anchor: int = _DEFAULT_FLANK_ANCHOR,
) -> list[Cluster]:
    """Pileup+cluster wrapper that returns ``[]`` instead of raising on empty loci."""
    counts = PileupCounts() if trace is not None else None
    try:
        obs = pileup_locus(
            bam,
            system.chromosome,
            system.ref_start - 1,
            system.ref_end,
            min_mapq=min_mapq,
            flank_anchor=flank_anchor,
            reference_fasta=reference_fasta,
            counts=counts,
        )
    except Exception:
        return []
    if trace is not None:
        trace.counts = counts
    if not obs:
        return []

    motifs = [m for m in system.motif.split(",") if m]
    bin_sizes: dict[int, int] | None = {} if trace is not None else None
    core_binned: dict[int, int] | None = {} if trace is not None else None
    clusters = cluster_observations(
        obs,
        identity_threshold=identity_threshold,
        # A per-marker override in the panel wins over the run-wide default.
        len_tolerance_bp=system.ont_len_tolerance or len_tolerance_bp,
        motifs=motifs,
        strand=system.strand,
        bin_sizes=bin_sizes,
        core_binned=core_binned,
    )
    if trace is not None and bin_sizes is not None and core_binned is not None:
        trace.binned_on_core = bool(motifs)
        trace.bins = [
            BinTrace(key_bp=key, n_reads=n, from_core=core_binned.get(key, 0) == n)
            for key, n in bin_sizes.items()
        ]
        trace.consensus_backend = poa_backend_name()
    return clusters


def _flag_non_default_thresholds(results: list[MarkerResult], params: RunParameters) -> None:
    """Mark every marker when the run overrode a *measured* default.

    Marker-level rather than run-level because that is where flags already live
    and where every consumer already looks — the audit census, the XLSX QC
    sheet, the HTML row tint. A run-level-only note would be the one thing a
    reviewer scanning per-locus rows never sees.

    Only ``derived`` provenance marks the run. Tuning a chosen threshold is
    ordinary work; overriding one computed from measured data means the profile
    is not comparable with a default run, and nobody remembers that six months
    later.
    """
    derived = params.derived_overrides()
    if not derived:
        return
    detail = "; ".join(f"{s.name}={params[s.name]} (default {s.default})" for s in derived)
    for result in results:
        if any(f.code == FlagCode.NON_DEFAULT_THRESHOLD for f in result.flags):
            continue
        result.flags.append(
            Flag.of(
                FlagCode.NON_DEFAULT_THRESHOLD,
                f"Run overrode {len(derived)} threshold(s) whose default was derived "
                f"from measured data: {detail}. This profile is not comparable with "
                "a default run.",
            )
        )


def _report_naming_fallbacks(system: System, alleles: list[Allele], called: list[Allele]) -> None:
    """Log every allele that fell back to the legacy CE, and how much it matters.

    Reported here rather than inside :mod:`frontstr.interp.naming` because the
    severity depends entirely on context the namer does not have. A one-read
    noise cluster failing the anchor identity guard is the guard doing its job —
    unpolished ONT reads carry indels that mangle a 30 bp flank — and belongs at
    debug. The same failure on a **called** allele means a reported number came
    from arithmetic known to be wrong at six markers, and has to be a warning.
    """
    called_ids = {id(a) for a in called}
    for a in alleles:
        if not a.strnaming_status or a.strnaming_status == NameStatus.OK.value:
            continue
        if a.strnaming_status == NameStatus.NO_RANGE.value:
            continue  # expected and permanent for DYS393/AMEL; not news
        event = "naming.fallback_on_called_allele" if id(a) in called_ids else "naming.fallback"
        emit = log.warning if id(a) in called_ids else log.debug
        emit(
            event,
            marker=system.name,
            reason=a.strnaming_status,
            n_reads=a.n_reads_total,
            consensus_bp=a.length_bp,
            consensus_method=a.consensus_method,
            number_method=a.number_method,
        )


def _allele_from_cluster(
    idx: int,
    c: Cluster,
    system: System,
    ref_length_bp: int,
    namer: StrNamer | None = None,
) -> Allele:
    """Build an unclassified :class:`Allele` from one :class:`Cluster`.

    The legacy CE is computed unconditionally: it is what the allele falls back
    to for markers STRNaming defines no range for (DYS393, AMEL) and whenever
    the reporting range cannot be located in the consensus.
    """
    consensus = c.consensus
    is_deletion = len(consensus) == 0
    isfg = compress_isfg(consensus, motif=system.motif, strand=system.strand) if consensus else ""
    if system.period == -1:
        raw_ce = ce_from_brackets(isfg) if isfg else None
        ce = (raw_ce - system.corr_value) if raw_ce is not None else None
    else:
        ce = ce_from_length(len(consensus), system.period, system.corr_value)
    bp_diff = len(consensus) - ref_length_bp if ref_length_bp is not None else 0
    allele_num, allele_src = compute_allele_numeric(len(consensus), system, ref_length_bp)

    named = namer.name(system.name, consensus) if namer is not None else None
    blocks = c.phase_sets
    return Allele(
        cluster_index=idx,
        consensus=consensus,
        strnaming_name=named.name if named else "",
        strnaming_ce=named.ce if named and named.ok else None,
        strnaming_status=named.status.value if named else NameStatus.NO_RANGE.value,
        length_bp=len(consensus),
        n_reads_total=c.n_reads,
        n_reads_hp1=c.n_hp1,
        n_reads_hp2=c.n_hp2,
        n_reads_hp_none=c.n_hp_none,
        phase_set=blocks.most_common(1)[0][0] if len(blocks) == 1 else None,
        n_phase_sets=len(blocks),
        n_forward=c.n_forward,
        n_reverse=c.n_reverse,
        mean_qual=c.mean_qual,
        ce=ce,
        isfg=isfg,
        bp_diff=bp_diff,
        is_deletion=is_deletion,
        consensus_method=c.consensus_method.value,
        allele_numeric=allele_num,
        allele_numeric_source=allele_src,
    )
