"""Cross-check between the Evidence layer (clusters) and LongTR (VCF).

We do *not* trust LongTR for forensic reporting — the counts come from our
evidence layer. But when LongTR is confident (``Q >= 0.9``) and its called
set of bp differences disagrees with the alleles we promoted, we raise a
**discordance flag**. Discordance never blocks a result; it surfaces in the
report as a yellow flag for analyst review.

Implements plan-longtr-improved.md §6.5.
"""

from __future__ import annotations

from frontstr.caller.vcf import LongTRResult
from frontstr.interp.models import Allele, AlleleStatus, Flag, FlagCode, MarkerResult

CONFIDENT_LONGTR_Q = 0.9


def cross_check(result: MarkerResult, longtr: LongTRResult | None) -> None:
    """Annotate ``result`` and its alleles with LongTR concordance flags.

    Mutates ``result`` and its alleles in place.
    """
    result.longtr_result = longtr
    if longtr is None:
        return

    inexact_sequences = {
        a.sequence for a in longtr.alleles if a.inexact and not a.is_deletion
    }
    for evidence_allele in result.alleles:
        for la in longtr.alleles:
            if _matches_longtr(evidence_allele, la):
                evidence_allele.longtr_match = True
                evidence_allele.longtr_bp_diff = la.bp_diff
                if la.inexact:
                    evidence_allele.longtr_inexact = True
                break

    if not longtr.samples:
        return
    sample = next(iter(longtr.samples.values()))
    if sample.gt_indices is None:
        return
    if sample.posterior is None or sample.posterior < CONFIDENT_LONGTR_Q:
        return

    longtr_bp_set = sorted(
        {longtr.alleles[i].bp_diff for i in sample.gt_indices if 0 <= i < len(longtr.alleles)}
    )
    evidence_bp_set = sorted(
        {
            a.bp_diff
            for a in result.alleles_called
            if a.status in (AlleleStatus.ALLELE, AlleleStatus.INEXACT_ALLELE)
        }
    )
    if longtr_bp_set and evidence_bp_set and longtr_bp_set != evidence_bp_set:
        result.discordant = True
        result.flags.append(
            Flag.of(
                FlagCode.LONGTR_DISCORDANT,
                f"LongTR called bp={longtr_bp_set} but evidence layer called "
                f"bp={evidence_bp_set}",
            )
        )

    # Propagate INEXACT_ALLELE upgrade for clusters that match an inexact LongTR ALT.
    for a in result.alleles:
        if a.status == AlleleStatus.ALLELE and a.consensus in inexact_sequences:
            a.status = AlleleStatus.INEXACT_ALLELE


def _matches_longtr(evidence: Allele, longtr_allele: object) -> bool:
    """Match an evidence allele to a LongTR allele.

    Exact sequence match preferred; falls back to ``bp_diff`` equality when
    sequences disagree (LongTR may have introduced ``INEXACT_ALLELE`` POA
    consensus that differs from ours).
    """
    seq = getattr(longtr_allele, "sequence", "")
    bp = getattr(longtr_allele, "bp_diff", None)
    is_del = getattr(longtr_allele, "is_deletion", False)
    if evidence.is_deletion and is_del:
        return True
    if evidence.consensus == seq:
        return True
    return bp is not None and bp == evidence.bp_diff
