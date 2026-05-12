"""Tri-allelic and mixture call rule.

Implements the call profile from plan-longtr-improved.md §6.4.

Outcomes:

- 0 candidates → ``no_data``
- 1 candidate, or 2nd allele below ``min_phr`` of the 1st → ``homozygous``
- 2 candidates (or 3rd below calling threshold) → ``heterozygous``
- 3+ candidates with ``allow_triallelic``:
    - all three balanced (PHR ≥ ``tri_balanced_thr``) → ``triallelic_type_II``
    - top two balanced, 3rd between 0.3 and ``tri_balanced_thr`` → ``triallelic_type_I``
    - any other 3-peak pattern → ``triallelic_review``
- 3+ candidates without ``allow_triallelic`` → ``two_called_three_present_review``
  (and ``tri_type = mixture_suspected``)

Type-II thresholds: default ``tri_balanced_thr = 0.6``; TPOX panels typically
configure 0.5 to capture balanced TPOX type-II patterns.
"""

from __future__ import annotations

from frontstr.interp.models import Allele, AlleleStatus, CallRule, TriType
from frontstr.panel.models import System

DEFAULT_TRI_BALANCED_THR = 0.6
DEFAULT_MIN_PHR_FOR_HET = 0.4  # 2nd allele must be >= 40% of 1st to call het
DEFAULT_TRI_TYPE_I_LOW = 0.3


def call_profile(
    alleles: list[Allele],
    system: System,
    *,
    calling_thresh: float = 0.10,
    min_phr_for_het: float = DEFAULT_MIN_PHR_FOR_HET,
    tri_type_i_low: float = DEFAULT_TRI_TYPE_I_LOW,
) -> tuple[list[Allele], CallRule, TriType]:
    """Decide which alleles to report and the overall call rule.

    Args:
        alleles: All classified candidates for the marker (any status).
        system: Marker definition; ``allow_triallelic`` and
            ``tri_balanced_thr`` are honoured.
        calling_thresh: Same value used by :mod:`classify`; used to decide
            whether a 3rd candidate is real or below the calling floor.
        min_phr_for_het: A 2nd allele below this fraction of the 1st is
            collapsed into homozygous.
        tri_type_i_low: Lower bound for the 3rd allele in a Type-I pattern.

    Returns:
        ``(alleles_called, call_rule, tri_type)``.
    """
    cands = sorted(
        (
            a
            for a in alleles
            if a.status in (AlleleStatus.ALLELE, AlleleStatus.INEXACT_ALLELE)
            and not a.is_deletion
        ),
        key=lambda a: a.n_reads_total,
        reverse=True,
    )
    if not cands:
        return [], CallRule.NO_DATA, TriType.NONE

    top = cands[0]
    tri_thr = system.tri_balanced_thr or DEFAULT_TRI_BALANCED_THR

    if len(cands) == 1:
        return [top], CallRule.HOMOZYGOUS, TriType.NONE

    a1, a2 = cands[0], cands[1]
    phr_12 = a2.n_reads_total / max(a1.n_reads_total, 1)
    if phr_12 < min_phr_for_het:
        return [a1], CallRule.HOMOZYGOUS, TriType.NONE

    if len(cands) == 2:
        return [a1, a2], CallRule.HETEROZYGOUS, TriType.NONE

    a3 = cands[2]
    phr_13 = a3.n_reads_total / max(a1.n_reads_total, 1)
    if phr_13 < calling_thresh:
        return [a1, a2], CallRule.HETEROZYGOUS, TriType.NONE

    if not system.allow_triallelic:
        return (
            [a1, a2],
            CallRule.TWO_CALLED_THREE_PRESENT,
            TriType.MIXTURE_SUSPECTED,
        )

    phr_23 = a3.n_reads_total / max(a2.n_reads_total, 1)
    if phr_12 >= tri_thr and phr_23 >= tri_thr:
        return [a1, a2, a3], CallRule.TRIALLELIC_TYPE_II, TriType.TYPE_II_BALANCED
    if phr_12 >= tri_thr and tri_type_i_low <= phr_23 < tri_thr:
        return [a1, a2, a3], CallRule.TRIALLELIC_TYPE_I, TriType.TYPE_I_UNBALANCED
    return [a1, a2, a3], CallRule.TRIALLELIC_REVIEW, TriType.MIXTURE_SUSPECTED
