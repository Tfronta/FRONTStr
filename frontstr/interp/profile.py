"""Orchestrator: clusters + LongTR → :class:`MarkerResult`.

This is the public entry point of the Interpretation layer. Higher-level
code (CLI, report, exports) calls :func:`interpret_marker` (single locus)
or :func:`interpret_run` (full panel).
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from frontstr.caller.vcf import LongTRResult
from frontstr.evidence.cluster import Cluster, cluster_observations
from frontstr.evidence.pileup import pileup_locus
from frontstr.interp.amel import interpret_amel
from frontstr.interp.classify import classify_allele
from frontstr.interp.concordance import cross_check
from frontstr.interp.isfg import ce_from_brackets, ce_from_length, compress_isfg
from frontstr.interp.models import Allele, MarkerResult
from frontstr.interp.stutter import build_expected_stutter
from frontstr.interp.triallelic import call_profile
from frontstr.panel.models import Panel, System

DEFAULT_ANALYTICAL_THRESH = 0.02
DEFAULT_CALLING_THRESH = 0.10
DEFAULT_PARENT_FRACTION = 0.20


def interpret_marker(
    *,
    system: System,
    clusters: list[Cluster],
    longtr: LongTRResult | None = None,
    analytical_thresh: float = DEFAULT_ANALYTICAL_THRESH,
    calling_thresh: float = DEFAULT_CALLING_THRESH,
    parent_fraction: float = DEFAULT_PARENT_FRACTION,
    ref_length_bp: int | None = None,
) -> MarkerResult:
    """Interpret one marker's evidence + LongTR call into a :class:`MarkerResult`.

    Args:
        system: Marker definition.
        clusters: Output of :func:`frontstr.evidence.cluster.cluster_observations`.
        longtr: Optional LongTR record for cross-checking.
        analytical_thresh: Below this fraction → noise.
        calling_thresh: Below this fraction (but >= analytical) → artefact.
        parent_fraction: Clusters with fraction ≥ this become parents for
            stutter-expectation calculation (default 0.20).
        ref_length_bp: Reference allele length used to compute ``bp_diff``.
            If ``None`` and ``longtr`` is provided, falls back to ``len(longtr.alleles[0].sequence)``.
    """
    total_reads = sum(c.n_reads for c in clusters)
    if ref_length_bp is None and longtr and longtr.alleles:
        ref_length_bp = len(longtr.alleles[0].sequence)

    alleles = [
        _allele_from_cluster(idx, c, system, ref_length_bp)
        for idx, c in enumerate(clusters)
    ]

    parents = [
        c for c, a in zip(clusters, alleles, strict=True)
        if a.fraction(total_reads) >= parent_fraction and not a.is_deletion
    ]
    expected = build_expected_stutter(parents, system)

    inexact_seqs = frozenset(
        a.sequence for a in (longtr.alleles if longtr else []) if a.inexact
    )

    for a in alleles:
        a.expected_stutter = expected.get(a.consensus, 0.0)
        a.status = classify_allele(
            a,
            total_reads=total_reads,
            expected_stutter=expected,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
            longtr_inexact_seqs=inexact_seqs,
        )

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
    cross_check(result, longtr)
    return result


def interpret_run(
    *,
    bam: Path,
    panel: Panel,
    longtr_results: dict[str, LongTRResult] | None = None,
    min_mapq: int = 20,
    identity_threshold: float = 0.97,
    len_tolerance_bp: int = 0,
    analytical_thresh: float = DEFAULT_ANALYTICAL_THRESH,
    calling_thresh: float = DEFAULT_CALLING_THRESH,
    reference_fasta: Path | None = None,
) -> list[MarkerResult]:
    """End-to-end: for each marker in ``panel``, pileup → cluster → interpret.

    Args:
        bam: Indexed sample BAM or CRAM.
        panel: Panel definition.
        longtr_results: Optional dict ``{marker_name: LongTRResult}`` for
            concordance checks. Markers without an entry are interpreted
            from evidence alone.
        reference_fasta: Reference FASTA path; required when ``bam`` is a CRAM.

    Returns:
        One :class:`MarkerResult` per marker in panel order. Empty pileups
        produce a ``NO_DATA`` result rather than failing.
    """
    longtr_results = longtr_results or {}
    out: list[MarkerResult] = []
    for system in panel.systems:
        if system.marker_type == "amel":
            out.append(interpret_amel(system, bam, min_mapq=min_mapq,
                                      reference_fasta=reference_fasta))
            continue
        clusters = _safe_pileup_and_cluster(
            bam=bam, system=system, min_mapq=min_mapq,
            identity_threshold=identity_threshold, len_tolerance_bp=len_tolerance_bp,
            reference_fasta=reference_fasta,
        )
        out.append(
            interpret_marker(
                system=system,
                clusters=clusters,
                longtr=longtr_results.get(system.name),
                analytical_thresh=analytical_thresh,
                calling_thresh=calling_thresh,
            )
        )
    return out


def _safe_pileup_and_cluster(
    *, bam: Path, system: System, min_mapq: int,
    identity_threshold: float, len_tolerance_bp: int,
    reference_fasta: Path | None = None,
) -> list[Cluster]:
    """Pileup+cluster wrapper that returns ``[]`` instead of raising on empty loci."""
    try:
        obs = pileup_locus(
            bam, system.chromosome, system.ref_start - 1, system.ref_end,
            min_mapq=min_mapq, reference_fasta=reference_fasta,
        )
    except Exception:
        return []
    if not obs:
        return []
    return cluster_observations(
        obs,
        identity_threshold=identity_threshold,
        len_tolerance_bp=len_tolerance_bp,
    )


def _allele_from_cluster(
    idx: int, c: Cluster, system: System, ref_length_bp: int | None,
) -> Allele:
    """Build an unclassified :class:`Allele` from one :class:`Cluster`."""
    consensus = c.consensus
    is_deletion = len(consensus) == 0
    isfg = compress_isfg(consensus, motif=system.motif, strand=system.strand) if consensus else ""
    if system.period == -1:
        ce = ce_from_brackets(isfg) if isfg else None
    else:
        ce = ce_from_length(len(consensus), system.period, system.corr_value)
    bp_diff = (
        len(consensus) - ref_length_bp
        if ref_length_bp is not None
        else 0
    )
    return Allele(
        cluster_index=idx,
        consensus=consensus,
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
    )


def index_longtr_results(results: Iterable[LongTRResult]) -> dict[str, LongTRResult]:
    """Convenience: turn a list of LongTRResult into a marker-name keyed map."""
    return {r.marker_name: r for r in results}
