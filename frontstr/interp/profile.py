"""Orchestrator: evidence clusters → :class:`MarkerResult`.

This is the public entry point of the Interpretation layer. Higher-level
code (CLI, report, exports) calls :func:`interpret_marker` (single locus)
or :func:`interpret_run` (full panel).
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from frontstr.evidence.cluster import Cluster, cluster_observations
from frontstr.evidence.pileup import pileup_locus
from frontstr.interp.allele_numeric import compute_allele_numeric, resolve_ref_anchor_bp
from frontstr.interp.amel import interpret_amel
from frontstr.interp.catalog import annotate_alleles
from frontstr.interp.classify import classify_allele
from frontstr.interp.flags import derive_marker_flags
from frontstr.interp.haplotype import suppress_hp_phantoms
from frontstr.interp.isfg import ce_from_brackets, ce_from_length, compress_isfg
from frontstr.interp.models import Allele, MarkerResult
from frontstr.interp.naming import NameStatus, StrNamer, default_namer
from frontstr.interp.qc import QcThresholds, derive_run_qc_flags
from frontstr.interp.stutter import build_expected_stutter
from frontstr.interp.triallelic import call_profile
from frontstr.log import get_logger
from frontstr.panel.catalog import AlleleCatalog
from frontstr.panel.models import Panel, System

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
    suppress_hp_phantoms(alleles, system)

    # Enrich ISFG / CE / iso-allele suffix from the curated catalog (no-op if None).
    annotate_alleles(alleles, system, catalog)

    alleles_called, call_rule, tri_type = call_profile(
        alleles, system, calling_thresh=calling_thresh
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
) -> list[MarkerResult]:
    """End-to-end: for each marker in ``panel``, pileup → cluster → interpret.

    Args:
        bam: Indexed sample BAM or CRAM.
        panel: Panel definition.
        reference_fasta: Reference FASTA path; required when ``bam`` is a CRAM.
        qc_thresholds: Laboratory QC policy for the run-level flags (coverage,
            strand bias). Defaults apply when omitted.

    Returns:
        One :class:`MarkerResult` per marker in panel order, each carrying its
        marker-level *and* run-level QC flags. Empty pileups produce a
        ``NO_DATA`` result rather than failing.
    """
    log = get_logger(__name__)
    out: list[MarkerResult] = []
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
            out.append(
                interpret_amel(system, bam, min_mapq=min_mapq, reference_fasta=reference_fasta)
            )
            continue
        clusters = _safe_pileup_and_cluster(
            bam=bam,
            system=system,
            min_mapq=min_mapq,
            identity_threshold=identity_threshold,
            len_tolerance_bp=len_tolerance_bp,
            reference_fasta=reference_fasta,
        )
        result = interpret_marker(
            system=system,
            clusters=clusters,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
            catalog=catalog,
            namer=namer,
        )
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
    log.info(
        "run.complete",
        n_markers=len(out),
        n_called=sum(1 for r in out if r.alleles_called),
        flags=dict(sorted(Counter(f.code.value for r in out for f in r.flags).items())),
        low_coverage_reads=thresholds.low_coverage_reads,
    )
    return out


def _safe_pileup_and_cluster(
    *,
    bam: Path,
    system: System,
    min_mapq: int,
    identity_threshold: float,
    len_tolerance_bp: int,
    reference_fasta: Path | None = None,
) -> list[Cluster]:
    """Pileup+cluster wrapper that returns ``[]`` instead of raising on empty loci."""
    try:
        obs = pileup_locus(
            bam,
            system.chromosome,
            system.ref_start - 1,
            system.ref_end,
            min_mapq=min_mapq,
            reference_fasta=reference_fasta,
        )
    except Exception:
        return []
    if not obs:
        return []
    return cluster_observations(
        obs,
        identity_threshold=identity_threshold,
        # A per-marker override in the panel wins over the run-wide default.
        len_tolerance_bp=system.ont_len_tolerance or len_tolerance_bp,
        motifs=[m for m in system.motif.split(",") if m],
        strand=system.strand,
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
