"""Tests for structured :class:`Flag` derivation."""

from __future__ import annotations

from frontstr.interp.flags import derive_marker_flags
from frontstr.interp.models import (
    CallRule,
    Flag,
    FlagCode,
    FlagSeverity,
    MarkerResult,
    TriType,
)
from frontstr.panel.models import System


def _system(name: str = "TPOX") -> System:
    return System(
        name=name,
        chromosome="chr2",
        ref_start=1,
        ref_end=50,
        motif="AATG",
        period=4,
    )


def _result(tri: TriType, rule: CallRule) -> MarkerResult:
    return MarkerResult(
        marker_name="TPOX",
        system=_system(),
        alleles=[],
        alleles_called=[],
        call_rule=rule,
        tri_type=tri,
        total_reads=80,
    )


def test_flag_of_defaults_severity_from_code() -> None:
    f = Flag.of(FlagCode.LONGTR_DISCORDANT, "x")
    assert f.severity == FlagSeverity.WARN
    assert f.code == FlagCode.LONGTR_DISCORDANT
    # Explicit severity overrides the default.
    assert Flag.of(FlagCode.TRIALLELIC, "y", FlagSeverity.ERROR).severity == FlagSeverity.ERROR


def test_derive_mixture_flag() -> None:
    r = _result(TriType.MIXTURE_SUSPECTED, CallRule.TWO_CALLED_THREE_PRESENT)
    derive_marker_flags(r)
    codes = {f.code for f in r.flags}
    assert FlagCode.MIXTURE_SUSPECTED in codes
    assert FlagCode.TRIALLELIC not in codes


def test_derive_triallelic_flag() -> None:
    r = _result(TriType.TYPE_II_BALANCED, CallRule.TRIALLELIC_TYPE_II)
    derive_marker_flags(r)
    assert {f.code for f in r.flags} == {FlagCode.TRIALLELIC}


def test_derive_is_idempotent_and_non_duplicating() -> None:
    r = _result(TriType.MIXTURE_SUSPECTED, CallRule.TWO_CALLED_THREE_PRESENT)
    derive_marker_flags(r)
    derive_marker_flags(r)
    assert sum(f.code == FlagCode.MIXTURE_SUSPECTED for f in r.flags) == 1


def test_no_flags_for_clean_het() -> None:
    r = _result(TriType.NONE, CallRule.HETEROZYGOUS)
    derive_marker_flags(r)
    assert r.flags == []


def _bare_allele(**kw: object):
    from frontstr.interp.models import Allele, AlleleStatus

    base: dict[str, object] = {
        "cluster_index": 0,
        "consensus": "",
        "length_bp": 0,
        "n_reads_total": 10,
        "n_reads_hp1": 5,
        "n_reads_hp2": 5,
        "n_reads_hp_none": 0,
        "n_forward": 5,
        "n_reverse": 5,
        "mean_qual": 30.0,
        "ce": None,
        "isfg": "",
        "bp_diff": 0,
        "is_deletion": False,
        "status": AlleleStatus.ALLELE,
    }
    base.update(kw)
    return Allele(**base)  # type: ignore[arg-type]


def test_canonical_number_period_ce() -> None:
    a = _bare_allele(ce=9.0, allele_numeric=9.0, allele_numeric_source="period_ce")
    assert (a.number, a.number_method) == (9.0, "period_ce")


def test_canonical_number_bracket_count_for_compound() -> None:
    a = _bare_allele(ce=13.0, allele_numeric=-3.0, allele_numeric_source="delta_only")
    assert (a.number, a.number_method) == (13.0, "bracket_count")


def test_canonical_number_delta_only_when_no_brackets() -> None:
    a = _bare_allele(ce=None, allele_numeric=-3.0, allele_numeric_source="delta_only")
    assert (a.number, a.number_method) == (-3.0, "delta")


def test_canonical_number_none_for_deletion() -> None:
    a = _bare_allele(
        ce=None,
        length_bp=0,
        is_deletion=True,
        allele_numeric=None,
        allele_numeric_source="deletion",
    )
    assert (a.number, a.number_method) == (None, "none")


def test_isoalleles_detected_by_same_number_diff_sequence() -> None:
    """Two called alleles, same number, different ISFG → ISOALLELE flag + marks."""
    a1 = _bare_allele(
        consensus="TCTA" * 15,
        length_bp=60,
        ce=15.0,
        isfg="[TCTA]15",
        allele_numeric=15.0,
        allele_numeric_source="period_ce",
    )
    a2 = _bare_allele(
        cluster_index=1,
        consensus="TCTG" + "TCTA" * 14,
        length_bp=60,
        ce=15.0,
        isfg="TCTG [TCTA]14",
        allele_numeric=15.0,
        allele_numeric_source="period_ce",
    )
    r = MarkerResult(
        marker_name="D3S1358",
        system=_system("D3S1358"),
        alleles=[a1, a2],
        alleles_called=[a1, a2],
        call_rule=CallRule.HETEROZYGOUS,
        tri_type=TriType.NONE,
        total_reads=20,
    )
    derive_marker_flags(r)
    assert any(f.code == FlagCode.ISOALLELE for f in r.flags)
    assert a1.iso.is_isoallele and a2.iso.is_isoallele


def test_no_isoallele_flag_for_distinct_numbers() -> None:
    a1 = _bare_allele(
        consensus="TCTA" * 12,
        length_bp=48,
        ce=12.0,
        isfg="[TCTA]12",
        allele_numeric=12.0,
        allele_numeric_source="period_ce",
    )
    a2 = _bare_allele(
        cluster_index=1,
        consensus="TCTA" * 14,
        length_bp=56,
        ce=14.0,
        isfg="[TCTA]14",
        allele_numeric=14.0,
        allele_numeric_source="period_ce",
    )
    r = MarkerResult(
        marker_name="CSF1PO",
        system=_system("CSF1PO"),
        alleles=[a1, a2],
        alleles_called=[a1, a2],
        call_rule=CallRule.HETEROZYGOUS,
        tri_type=TriType.NONE,
        total_reads=20,
    )
    derive_marker_flags(r)
    assert not any(f.code == FlagCode.ISOALLELE for f in r.flags)


def test_amel_like_none_numbers_not_isoallele() -> None:
    """Two number-less alleles (AMEL X/Y) must NOT be grouped as iso-alleles."""
    x = _bare_allele(
        consensus="", isfg="X", ce=None, allele_numeric=None, allele_numeric_source="unavailable"
    )
    y = _bare_allele(
        cluster_index=1,
        consensus="",
        isfg="Y",
        ce=None,
        allele_numeric=None,
        allele_numeric_source="unavailable",
    )
    r = MarkerResult(
        marker_name="AMEL",
        system=_system("AMEL"),
        alleles=[x, y],
        alleles_called=[x, y],
        call_rule=CallRule.HETEROZYGOUS,
        tri_type=TriType.NONE,
        total_reads=20,
    )
    derive_marker_flags(r)
    assert not any(f.code == FlagCode.ISOALLELE for f in r.flags)


# ---------------------------------------------------------------------------
# CONSENSUS_FALLBACK — an unpolished consensus must never pass silently
# ---------------------------------------------------------------------------


def _result_with(alleles: list) -> MarkerResult:
    return MarkerResult(
        marker_name="TH01",
        system=_system("TH01"),
        alleles=alleles,
        alleles_called=alleles,
        call_rule=CallRule.HETEROZYGOUS,
        tri_type=TriType.NONE,
        total_reads=40,
    )


def test_unpolished_consensus_raises_a_warn_flag() -> None:
    r = _result_with([_bare_allele(consensus_method="mode", isfg="[AATG]7")])
    derive_marker_flags(r)
    flag = next(f for f in r.flags if f.code == FlagCode.CONSENSUS_FALLBACK)
    assert flag.severity == FlagSeverity.WARN
    assert "pyspoa" in flag.message


def test_poa_consensus_raises_no_fallback_flag() -> None:
    r = _result_with([_bare_allele(consensus_method="poa_spoa", isfg="[AATG]7")])
    derive_marker_flags(r)
    assert FlagCode.CONSENSUS_FALLBACK not in {f.code for f in r.flags}


def test_single_read_consensus_raises_no_fallback_flag() -> None:
    """A 1-read cluster is unpolishable, not mis-configured — different signal."""
    r = _result_with([_bare_allele(consensus_method="single", isfg="[AATG]7")])
    derive_marker_flags(r)
    assert FlagCode.CONSENSUS_FALLBACK not in {f.code for f in r.flags}


def test_hp_phantom_collapse_is_flagged_for_audit() -> None:
    from frontstr.interp.models import AlleleStatus

    owner = _bare_allele(consensus_method="poa_spoa")
    phantom = _bare_allele(
        cluster_index=1,
        length_bp=252,
        n_reads_total=5,
        status=AlleleStatus.HP_PHANTOM,
        consensus_method="poa_spoa",
    )
    r = MarkerResult(
        marker_name="D8S1179",
        system=_system("D8S1179"),
        alleles=[owner, phantom],
        alleles_called=[owner],
        call_rule=CallRule.HETEROZYGOUS,
        tri_type=TriType.NONE,
        total_reads=30,
    )
    derive_marker_flags(r)
    flag = next(f for f in r.flags if f.code == FlagCode.HP_PHANTOM_COLLAPSED)
    assert flag.severity == FlagSeverity.INFO
    assert "252 bp / 5 reads" in flag.message
