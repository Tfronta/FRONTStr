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
        name=name, chromosome="chr2", ref_start=1, ref_end=50, motif="AATG", period=4,
    )


def _result(tri: TriType, rule: CallRule) -> MarkerResult:
    return MarkerResult(
        marker_name="TPOX", system=_system(), alleles=[], alleles_called=[],
        call_rule=rule, tri_type=tri, total_reads=80,
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
