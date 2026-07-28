"""Unit tests for :mod:`frontstr.interp.allele_numeric`."""

from __future__ import annotations

from frontstr.interp.allele_numeric import (
    compute_allele_numeric,
    resolve_ref_anchor_bp,
)
from frontstr.panel.models import System


def _tetramer(name: str = "TH01_X") -> System:
    return System(
        name=name,
        chromosome="chr11",
        ref_start=1,
        ref_end=48,
        motif="AATG",
        period=4,
        corr_value=0,
    )


def test_resolve_explicit_overrides_panel_span() -> None:
    s = _tetramer()
    assert resolve_ref_anchor_bp(s, explicit=99) == 99
    assert s.span() == 48


def test_period_ce_simple_locus() -> None:
    s = _tetramer()
    n, src = compute_allele_numeric(48, s, ref_anchor_bp=48)
    assert src == "period_ce"
    assert n == 12.0


def test_reference_offset_compound_marker() -> None:
    s = System(
        name="D3X",
        chromosome="chr3",
        ref_start=100,
        ref_end=163,
        motif="TCTA,TCTG",
        period=-1,
        corr_value=0,
        reference_ce=15.0,
        allele_bp_step=4,
    )
    anchor = s.span()
    assert anchor == 64
    n, src = compute_allele_numeric(68, s, ref_anchor_bp=anchor)
    assert src == "reference_offset"
    assert n == 16.0


def test_delta_only_without_reference_ce() -> None:
    s = System(
        name="D3Y",
        chromosome="chr3",
        ref_start=100,
        ref_end=163,
        motif="TCTA,TCTG",
        period=-1,
        corr_value=0,
        reference_ce=None,
        allele_bp_step=4,
    )
    anchor = s.span()
    n, src = compute_allele_numeric(68, s, ref_anchor_bp=anchor)
    assert src == "delta_only"
    assert n == 1.0
