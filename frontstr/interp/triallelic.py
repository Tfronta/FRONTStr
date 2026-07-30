"""Tri-allelic and mixture call rule.

Implements the call profile from plan-longtr-improved.md §6.4.

Outcomes:

- 0 candidates → ``no_data``
- 1 candidate, or 2nd allele below ``min_phr`` of the 1st → ``homozygous``,
  **unless phasing puts the two on opposite haplotypes** — see below
- 2 candidates (or 3rd below calling threshold, or 3rd below the absolute
  ``min_reads_third`` floor) → ``heterozygous``
- 3+ candidates with ``allow_triallelic``:
    - all three balanced (PHR ≥ ``tri_balanced_thr``) → ``triallelic_type_II``
    - top two balanced, 3rd between 0.3 and ``tri_balanced_thr`` → ``triallelic_type_I``
    - any other 3-peak pattern → ``triallelic_review``
- 3+ candidates without ``allow_triallelic`` → ``two_called_three_present_review``
  (and ``tri_type = mixture_suspected``)

Type-II thresholds: default ``tri_balanced_thr = 0.6``; TPOX panels typically
configure 0.5 to capture balanced TPOX type-II patterns.

The peak-height ratio is not the best evidence available
--------------------------------------------------------

``min_phr_for_het`` is inherited from capillary electrophoresis, where peak
height is all there is. On phased long reads it is not: a diploid locus carries
one allele per haplotype, so two candidates confidently on opposite haplotypes
are two alleles regardless of their read ratio. At ONT depths a 5-vs-17 split
(PHR 0.29) is ordinary sampling, and collapsing it to a homozygote is a
false exclusion waiting to happen.

So the PHR floor is now overridden by
:func:`frontstr.interp.haplotype.on_opposite_haplotypes`, and the rescued
allele carries ``hp_rescued`` so the marker can raise ``HP_RESCUED_HET``. This
is the same invariant :func:`~frontstr.interp.haplotype.suppress_hp_phantoms`
uses to *delete* alleles, applied in the opposite direction. Both are no-ops on
unphased BAMs, where the ratio is again the only evidence there is.
"""

from __future__ import annotations

from frontstr.interp.haplotype import on_opposite_haplotypes
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
    min_reads_third: int | None = None,
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
        min_reads_third: Absolute read-count floor a 3rd candidate must clear
            (in addition to ``calling_thresh``) before it can promote the locus
            to triallelic or raise ``mixture_suspected``. Defaults to
            ``system.min_reads_third``. Guards against ONT basecaller phantoms
            (known-bug #6).

    Returns:
        ``(alleles_called, call_rule, tri_type)``.
    """
    cands = sorted(
        (
            a
            for a in alleles
            if a.status in (AlleleStatus.ALLELE, AlleleStatus.INEXACT_ALLELE) and not a.is_deletion
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
        # The peak-height ratio is an indirect proxy for "these are two
        # alleles"; phasing is direct evidence of it. A diploid locus carries
        # one allele per haplotype, so two candidates confidently on opposite
        # haplotypes are two alleles no matter how unbalanced their coverage —
        # and at ONT depths an unlucky 5-vs-17 split is ordinary sampling, not
        # evidence of homozygosity. No-op on unphased BAMs, where the ratio
        # remains the only evidence available.
        if not on_opposite_haplotypes(a1, a2):
            return [a1], CallRule.HOMOZYGOUS, TriType.NONE
        a2.hp_rescued = True

    if len(cands) == 2:
        return [a1, a2], CallRule.HETEROZYGOUS, TriType.NONE

    a3 = cands[2]
    phr_13 = a3.n_reads_total / max(a1.n_reads_total, 1)
    read_floor = system.min_reads_third if min_reads_third is None else min_reads_third
    if phr_13 < calling_thresh or a3.n_reads_total < read_floor:
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
