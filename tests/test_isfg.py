"""Tests for ISFG bracketed nomenclature compression."""

from __future__ import annotations

from frontstr.interp.isfg import ce_from_length, compress_isfg, motif_repeat_summary


def test_single_motif_pure_run() -> None:
    assert compress_isfg("AGATAGATAGAT", motif="AGAT") == "[AGAT]3"


def test_single_motif_with_flank() -> None:
    out = compress_isfg("TAGATAGATAGATC", motif="AGAT")
    assert "[AGAT]3" in out


def test_compress_isfg_collapses_non_motif_bases() -> None:
    """Non-motif stretches must not be one space per nucleotide."""
    seq = "CCGAT" + "TCTA" * 3 + "TT"
    out = compress_isfg(seq, motif="TCTA,TCTG")
    assert " C " not in out and " T " not in out
    assert out == "CCGAT [TCTA]3 TT"


def test_multi_motif_d3s1358() -> None:
    seq = "TCTA" * 4 + "TCTG" * 1 + "TCTA" * 3
    out = compress_isfg(seq, motif="TCTA,TCTG")
    assert "[TCTA]4" in out
    assert "[TCTA]3" in out


def test_motif_repeat_summary_d3_like() -> None:
    seq = "TCTA" * 4 + "TCTG" + "TCTA" * 2
    s = motif_repeat_summary(seq, "TCTA,TCTG")
    assert "TCTA" in s and "TCTG" in s
    assert "TR" in s and "bp" in s


def test_single_repeat_unit_not_bracketed() -> None:
    """A single occurrence of the motif should not be bracketed."""
    out = compress_isfg("AGAT", motif="AGAT")
    assert out == "AGAT"


def test_empty_input() -> None:
    assert compress_isfg("", motif="AGAT") == ""


def test_no_motif_provided() -> None:
    assert compress_isfg("ACGTACGT", motif="") == "ACGTACGT"


def test_ce_from_length_basic() -> None:
    assert ce_from_length(48, period=4, corr_value=0) == 12.0


def test_ce_from_length_with_corr() -> None:
    assert ce_from_length(48, period=4, corr_value=4) == 11.0


def test_ce_from_length_microvariant() -> None:
    assert ce_from_length(49, period=4, corr_value=0) == 12.2


def test_ce_undefined_for_multi_motif() -> None:
    assert ce_from_length(60, period=-1, corr_value=0) is None
