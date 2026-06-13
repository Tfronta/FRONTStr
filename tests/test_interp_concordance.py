"""Tests for evidence-vs-LongTR concordance checks."""

from __future__ import annotations

from frontstr.caller.vcf import LongTRAlleleSpec, LongTRResult, LongTRSampleCall
from frontstr.interp.concordance import cross_check
from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    FlagCode,
    FlagSeverity,
    MarkerResult,
    TriType,
)
from frontstr.panel.models import System


def _system() -> System:
    return System(
        name="TH01", chromosome="chr11", ref_start=2_171_000, ref_end=2_171_050,
        motif="AATG", period=4,
    )


def _allele(idx: int, consensus: str, n_reads: int, bp_diff: int,
            status: AlleleStatus = AlleleStatus.ALLELE) -> Allele:
    return Allele(
        cluster_index=idx, consensus=consensus, length_bp=len(consensus),
        n_reads_total=n_reads, n_reads_hp1=0, n_reads_hp2=0, n_reads_hp_none=n_reads,
        n_forward=n_reads, n_reverse=0, mean_qual=30.0,
        ce=float(len(consensus) / 4), isfg=f"[AATG]{len(consensus) // 4}",
        bp_diff=bp_diff, is_deletion=False, status=status,
    )


def _make_result(alleles: list[Allele], called: list[Allele]) -> MarkerResult:
    return MarkerResult(
        marker_name="TH01", system=_system(),
        alleles=alleles, alleles_called=called,
        call_rule=CallRule.HETEROZYGOUS, tri_type=TriType.NONE,
        total_reads=sum(a.n_reads_total for a in alleles),
    )


def test_cross_check_none_longtr_noop() -> None:
    r = _make_result([], [])
    cross_check(r, None)
    assert r.longtr_result is None
    assert r.discordant is False


def test_cross_check_concordant_call() -> None:
    """Same bp set called by LongTR (high Q) and evidence → not discordant."""
    a1 = _allele(0, "AATG" * 8, 60, bp_diff=4)
    a2 = _allele(1, "AATG" * 7, 55, bp_diff=0)
    result = _make_result([a1, a2], [a1, a2])

    ref = "AATG" * 7
    alt = "AATG" * 8
    longtr = LongTRResult(
        marker_name="TH01", chrom="chr11", pos=2_171_000,
        motif="AATG", period="4",
        alleles=[
            LongTRAlleleSpec(sequence=ref, bp_diff=0, inexact=False, is_deletion=False),
            LongTRAlleleSpec(sequence=alt, bp_diff=4, inexact=False, is_deletion=False),
        ],
        samples={"S1": LongTRSampleCall(
            sample="S1", gt_indices=(0, 1), phased=False, posterior=0.99, depth=115,
            pdp_hp1=60, pdp_hp2=55,
        )},
    )
    cross_check(result, longtr)
    assert result.longtr_result is longtr
    assert result.discordant is False
    # Both evidence alleles should now have longtr_match=True
    assert a1.longtr_match is True
    assert a2.longtr_match is True


def test_cross_check_discordant_high_q() -> None:
    """LongTR confident but called a different bp set than evidence."""
    a1 = _allele(0, "AATG" * 8, 60, bp_diff=4)
    a2 = _allele(1, "AATG" * 7, 55, bp_diff=0)
    result = _make_result([a1, a2], [a1, a2])

    # LongTR called REF (0) + a 12-bp diff that's not in our evidence
    longtr = LongTRResult(
        marker_name="TH01", chrom="chr11", pos=2_171_000,
        motif="AATG", period="4",
        alleles=[
            LongTRAlleleSpec(sequence="AATG" * 7, bp_diff=0, inexact=False, is_deletion=False),
            LongTRAlleleSpec(sequence="AATG" * 10, bp_diff=12, inexact=False, is_deletion=False),
        ],
        samples={"S1": LongTRSampleCall(
            sample="S1", gt_indices=(0, 1), phased=False, posterior=0.95, depth=115,
            pdp_hp1=60, pdp_hp2=55,
        )},
    )
    cross_check(result, longtr)
    assert result.discordant is True
    assert any(f.code == FlagCode.LONGTR_DISCORDANT for f in result.flags)
    disc = next(f for f in result.flags if f.code == FlagCode.LONGTR_DISCORDANT)
    assert disc.severity == FlagSeverity.WARN
    assert "LongTR called bp" in disc.message


def test_cross_check_low_q_does_not_flag() -> None:
    """LongTR with Q < 0.9 must not raise the flag."""
    a1 = _allele(0, "AATG" * 8, 60, bp_diff=4)
    result = _make_result([a1], [a1])

    longtr = LongTRResult(
        marker_name="TH01", chrom="chr11", pos=2_171_000,
        motif="AATG", period="4",
        alleles=[
            LongTRAlleleSpec(sequence="AATG" * 7, bp_diff=0, inexact=False, is_deletion=False),
            LongTRAlleleSpec(sequence="AATG" * 10, bp_diff=12, inexact=False, is_deletion=False),
        ],
        samples={"S1": LongTRSampleCall(
            sample="S1", gt_indices=(0, 1), phased=False, posterior=0.5, depth=60,
        )},
    )
    cross_check(result, longtr)
    assert result.discordant is False


def test_cross_check_inexact_propagates() -> None:
    """An evidence allele matching a LongTR INEXACT_ALLELE gets upgraded."""
    a1 = _allele(0, "AATG" * 8, 60, bp_diff=4, status=AlleleStatus.ALLELE)
    result = _make_result([a1], [a1])

    longtr = LongTRResult(
        marker_name="TH01", chrom="chr11", pos=2_171_000,
        motif="AATG", period="4",
        alleles=[
            LongTRAlleleSpec(sequence="AATG" * 7, bp_diff=0, inexact=False, is_deletion=False),
            LongTRAlleleSpec(sequence="AATG" * 8, bp_diff=4, inexact=True, is_deletion=False),
        ],
        samples={"S1": LongTRSampleCall(
            sample="S1", gt_indices=(0, 1), phased=False, posterior=0.99, depth=60,
        )},
    )
    cross_check(result, longtr)
    assert a1.status == AlleleStatus.INEXACT_ALLELE
    assert a1.longtr_inexact is True
