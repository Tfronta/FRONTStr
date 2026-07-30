"""Tests for :mod:`frontstr.interp.qc` — the threshold-dependent QC flags.

These five codes were declared and never emitted for the whole life of the
project, so the tests care as much about *when a flag must not fire* as about
when it must: a QC layer that cries wolf on half the panel trains reviewers to
ignore it, which is worse than no layer at all.
"""

from __future__ import annotations

import pytest

from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    Flag,
    FlagCode,
    FlagSeverity,
    MarkerResult,
    TriType,
)
from frontstr.interp.qc import QcThresholds, derive_run_qc_flags, two_sided_binomial_p
from frontstr.panel.models import System


def _system(name: str = "TH01", **kw: object) -> System:
    base: dict[str, object] = {
        "name": name,
        "chromosome": "chr11",
        "ref_start": 1,
        "ref_end": 400,
        "motif": "AATG",
        "period": 4,
    }
    base.update(kw)
    return System(**base)  # type: ignore[arg-type]


def _allele(
    *,
    n_forward: int = 10,
    n_reverse: int = 10,
    status: AlleleStatus = AlleleStatus.ALLELE,
    ce: float | None = 9.0,
) -> Allele:
    n = n_forward + n_reverse
    return Allele(
        cluster_index=0,
        consensus="AATG" * 9,
        length_bp=36,
        n_reads_total=n,
        n_reads_hp1=0,
        n_reads_hp2=0,
        n_reads_hp_none=n,
        n_forward=n_forward,
        n_reverse=n_reverse,
        mean_qual=30.0,
        ce=ce,
        isfg="[AATG]9",
        bp_diff=0,
        is_deletion=False,
        status=status,
        allele_numeric=ce,
        allele_numeric_source="period_ce" if ce is not None else "",
    )


def _result(
    *,
    system: System | None = None,
    alleles: list[Allele] | None = None,
    call_rule: CallRule = CallRule.HETEROZYGOUS,
    total_reads: int = 40,
) -> MarkerResult:
    # Coverage is expressed through the alleles when none are given: a locus
    # cannot have 12 spanning reads and a 20-read allele, and low_coverage now
    # measures the reads supporting the call rather than the spanning total.
    called = alleles if alleles is not None else [_a(total_reads)]
    return MarkerResult(
        marker_name=(system or _system()).name,
        system=system or _system(),
        alleles=called,
        alleles_called=called,
        call_rule=call_rule,
        tri_type=TriType.NONE,
        total_reads=total_reads,
    )


def _codes(r: MarkerResult) -> set[FlagCode]:
    return {f.code for f in r.flags}


# ---------------------------------------------------------------------------
# two_sided_binomial_p
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("k", "n", "expected"),
    [
        (10, 20, 1.0),  # perfectly balanced
        (0, 10, 2 / 1024),  # both tails of an all-one-strand split
        (20, 20, 2 / 2**20),
    ],
)
def test_binomial_p_known_values(k: int, n: int, expected: float) -> None:
    assert two_sided_binomial_p(k, n) == pytest.approx(expected)


def test_binomial_p_is_symmetric() -> None:
    assert two_sided_binomial_p(3, 20) == two_sided_binomial_p(17, 20)


def test_binomial_p_of_no_reads_is_one() -> None:
    assert two_sided_binomial_p(0, 0) == 1.0


# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------


def test_no_call_is_a_dropout_not_low_coverage() -> None:
    """A locus that produced nothing is a different problem from a thin call."""
    r = _result(alleles=[], call_rule=CallRule.NO_DATA, total_reads=0)
    derive_run_qc_flags([r])
    assert _codes(r) == {FlagCode.DROPOUT}


def test_thin_call_is_low_coverage() -> None:
    r = _result(total_reads=12)
    derive_run_qc_flags([r])
    assert FlagCode.LOW_COVERAGE in _codes(r)
    assert FlagCode.DROPOUT not in _codes(r)


def test_adequate_coverage_raises_nothing() -> None:
    r = _result(total_reads=40)
    derive_run_qc_flags([r])
    assert r.flags == []


def test_low_coverage_floor_is_configurable() -> None:
    r = _result(total_reads=25)
    derive_run_qc_flags([r], QcThresholds(low_coverage_reads=30))
    assert FlagCode.LOW_COVERAGE in _codes(r)


def test_haploid_locus_gets_no_dropout_wording() -> None:
    """A Y marker has no second allele to lose; the generic wording would lie."""
    r = _result(system=_system("DYS391", category="y_chromosomal"), total_reads=12)
    derive_run_qc_flags([r])
    msg = next(f.message for f in r.flags if f.code == FlagCode.LOW_COVERAGE)
    assert "minor allele" not in msg
    assert "thin evidence" in msg


# ---------------------------------------------------------------------------
# Strand bias
# ---------------------------------------------------------------------------


def test_balanced_strands_raise_nothing() -> None:
    r = _result(alleles=[_allele(n_forward=15, n_reverse=15)])
    derive_run_qc_flags([r])
    assert FlagCode.STRAND_BIAS not in _codes(r)


def test_extreme_strand_skew_is_flagged() -> None:
    r = _result(alleles=[_allele(n_forward=18, n_reverse=1)])
    derive_run_qc_flags([r])
    assert FlagCode.STRAND_BIAS in _codes(r)


def test_small_alleles_are_not_strand_tested() -> None:
    """Below the minimum the test cannot reach the threshold even for 0/n.

    Testing anyway would produce a 'passed' that means nothing.
    """
    r = _result(alleles=[_allele(n_forward=6, n_reverse=0)], total_reads=40)
    derive_run_qc_flags([r])
    assert FlagCode.STRAND_BIAS not in _codes(r)


def test_strand_bias_message_names_the_allele_and_the_p_value() -> None:
    r = _result(alleles=[_allele(n_forward=20, n_reverse=1)])
    derive_run_qc_flags([r])
    msg = next(f.message for f in r.flags if f.code == FlagCode.STRAND_BIAS)
    assert "20+/1-" in msg
    assert "p=" in msg


# ---------------------------------------------------------------------------
# Inexact and kit nomenclature
# ---------------------------------------------------------------------------


def test_inexact_called_allele_is_flagged() -> None:
    r = _result(alleles=[_allele(status=AlleleStatus.INEXACT_ALLELE)])
    derive_run_qc_flags([r])
    assert FlagCode.INEXACT_ALLELE in _codes(r)


def test_kit_nomenclature_note_becomes_a_warning() -> None:
    system = _system("vWA", kit_nomenclature_note="Offsets run in opposite directions.")
    r = _result(system=system)
    derive_run_qc_flags([r])
    flag = next(f for f in r.flags if f.code == FlagCode.CE_NOMENCLATURE_OFFSET)
    assert flag.severity == FlagSeverity.WARN, "a silent offset can cause a false exclusion"
    assert "opposite directions" in flag.message


def test_marker_without_a_note_is_not_flagged() -> None:
    r = _result()
    derive_run_qc_flags([r])
    assert FlagCode.CE_NOMENCLATURE_OFFSET not in _codes(r)


def test_uncalled_marker_gets_no_nomenclature_warning() -> None:
    """Nothing was reported, so there is no number to misread."""
    system = _system("vWA", kit_nomenclature_note="…")
    r = _result(system=system, alleles=[], call_rule=CallRule.NO_DATA, total_reads=0)
    derive_run_qc_flags([r])
    assert FlagCode.CE_NOMENCLATURE_OFFSET not in _codes(r)


# ---------------------------------------------------------------------------
# Contract
# ---------------------------------------------------------------------------


def test_derivation_is_idempotent() -> None:
    r = _result(total_reads=5)
    derive_run_qc_flags([r])
    derive_run_qc_flags([r])
    assert sum(1 for f in r.flags if f.code == FlagCode.LOW_COVERAGE) == 1


def test_existing_flags_are_preserved() -> None:
    r = _result(total_reads=5)
    r.flags.append(Flag.of(FlagCode.MIXTURE_SUSPECTED, "from an earlier producer"))
    derive_run_qc_flags([r])
    assert FlagCode.MIXTURE_SUSPECTED in _codes(r)
    assert FlagCode.LOW_COVERAGE in _codes(r)


def test_returns_the_thresholds_actually_applied() -> None:
    """The audit record needs the effective policy, not the caller's intent."""
    applied = derive_run_qc_flags([_result()], QcThresholds(low_coverage_reads=17))
    assert applied.low_coverage_reads == 17
    assert derive_run_qc_flags([_result()]).low_coverage_reads == 20


# ---------------------------------------------------------------------------
# Allele balance
# ---------------------------------------------------------------------------


def _a(n_reads: int) -> Allele:
    """A called allele carrying exactly ``n_reads``."""
    return _allele(n_forward=n_reads // 2, n_reverse=n_reads - n_reads // 2)


def _het(*counts: int) -> MarkerResult:
    return _result(alleles=[_a(n) for n in counts], total_reads=sum(counts))


class TestAlleleBalance:
    """AB replaces the peak-height ratio in output; it must not contradict it."""

    def test_perfectly_balanced_het_is_half(self) -> None:
        assert _het(20, 20).allele_balance == 0.5

    def test_scale_runs_from_the_strongest_allele(self) -> None:
        """Stated convention: strongest over the pair, so AB is one-sided."""
        assert _het(17, 16).allele_balance == 0.515
        # The order of the called list must not change the answer.
        assert _het(16, 17).allele_balance == _het(17, 16).allele_balance

    def test_undefined_for_a_homozygote(self) -> None:
        assert _het(30).allele_balance is None

    def test_undefined_for_a_triallelic_locus(self) -> None:
        """One ratio across three alleles would be a fiction."""
        assert _het(20, 18, 17).allele_balance is None

    def test_undefined_when_no_allele_is_called(self) -> None:
        assert _result(alleles=[]).allele_balance is None

    def test_the_calling_floor_lands_outside_the_balanced_band(self) -> None:
        """min_phr_for_het = 0.4 is AB 0.714, so the band sits *inside* the
        callable range: it warns about a het rather than rejecting it."""
        assert QcThresholds().balanced_ab_max < 1 / (1 + 0.4)


class TestAlleleImbalanceFlag:
    def test_uneven_het_is_flagged(self) -> None:
        r = _het(20, 8)
        derive_run_qc_flags([r])
        assert FlagCode.ALLELE_IMBALANCE in _codes(r)

    def test_balanced_het_is_not_flagged(self) -> None:
        r = _het(17, 16)
        derive_run_qc_flags([r])
        assert FlagCode.ALLELE_IMBALANCE not in _codes(r)

    def test_a_phasing_rescue_suppresses_the_duplicate_warning(self) -> None:
        """HP_RESCUED_HET already says the read ratio was the problem.

        Two warnings for one phenomenon trains a reviewer to skim them.
        """
        r = _het(17, 5)
        r.alleles_called[1].hp_rescued = True
        assert r.allele_balance is not None
        assert r.allele_balance > QcThresholds().balanced_ab_max
        derive_run_qc_flags([r])
        assert FlagCode.ALLELE_IMBALANCE not in _codes(r)

    def test_homozygote_never_raises_it(self) -> None:
        r = _het(30)
        derive_run_qc_flags([r])
        assert FlagCode.ALLELE_IMBALANCE not in _codes(r)


class TestCoverageIsMeasuredOnTheCall:
    """The floor watches the evidence behind the genotype, not the pileup size.

    Reads that clustered into neither allele are not draws from the pair the
    binomial models, so counting them made the flag looser than derived — and
    it stayed silent on HG00263 D18S51, which is called on 11 reads out of 33
    spanning and misses the second allele Illumina sees.
    """

    def test_a_thinly_supported_call_is_flagged_despite_deep_pileup(self) -> None:
        r = _result(alleles=[_a(11)], total_reads=33)
        derive_run_qc_flags([r])
        assert FlagCode.LOW_COVERAGE in _codes(r)

    def test_a_well_supported_call_is_not_flagged(self) -> None:
        r = _result(alleles=[_a(15), _a(15)], total_reads=40)
        derive_run_qc_flags([r])
        assert FlagCode.LOW_COVERAGE not in _codes(r)

    def test_the_message_quotes_both_numbers(self) -> None:
        """A reviewer needs to see the gap, not just the smaller number."""
        r = _result(alleles=[_a(11)], total_reads=33)
        derive_run_qc_flags([r])
        msg = next(f.message for f in r.flags if f.code == FlagCode.LOW_COVERAGE)
        assert "11" in msg and "33" in msg

    def test_the_floor_sits_at_the_knee_of_the_risk_curve(self) -> None:
        """Guards the re-derivation: 20 is chosen because 17 -> 20 halves the
        dropout risk while 20 -> 25 barely moves it. If someone changes the
        default, the reasoning in QcThresholds has to be redone, not assumed.
        """
        from math import comb

        p_minor, calling, ratio = 0.4 / 1.4, 0.10, 0.795

        def worst_risk_from(floor: int) -> float:
            def risk(n: int) -> float:
                need = calling * (n / ratio)
                k_min = int(need) + (0 if float(need).is_integer() else 1)
                return sum(comb(n, k) * p_minor**k * (1 - p_minor) ** (n - k) for k in range(k_min))

            return max(risk(n) for n in range(floor, 81))

        assert QcThresholds().low_coverage_reads == 20
        assert worst_risk_from(17) > 0.09, "below the knee the risk climbs"
        assert worst_risk_from(20) < 0.06, "at the knee it is ~5.7%"
        assert worst_risk_from(25) > 0.04, "above it, more flags buy almost nothing"
