"""Tests for ISFG bracketed nomenclature compression."""

from __future__ import annotations

from frontstr.interp.isfg import (
    ce_from_brackets,
    ce_from_length,
    compress_isfg,
    motif_repeat_summary,
)


def test_single_motif_pure_run() -> None:
    assert compress_isfg("AGATAGATAGAT", motif="AGAT") == "[AGAT]3"


def test_single_motif_with_flank() -> None:
    out = compress_isfg("TAGATAGATAGATC", motif="AGAT")
    assert "[AGAT]3" in out


def test_compress_isfg_collapses_non_motif_bases() -> None:
    """Non-motif stretches must not be one space per nucleotide, and must be lowercase."""
    seq = "CCGAT" + "TCTA" * 3 + "TT"
    out = compress_isfg(seq, motif="TCTA,TCTG")
    assert " C " not in out and " T " not in out
    assert out == "ccgat [TCTA]3 tt"


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


def test_ce_from_length_microvariant_1bp() -> None:
    # 12 full repeats × 4 + 1 extra base  →  CE 12.1  (NOT 12.25 rounded)
    assert ce_from_length(49, period=4, corr_value=0) == 12.1


def test_ce_from_length_microvariant_3bp() -> None:
    # TH01 9.3: 9 full AATG (36 bp) + 3 extra bp (ATG) = 39 bp, corr_value=0
    assert ce_from_length(39, period=4, corr_value=0) == 9.3


def test_ce_from_length_microvariant_2bp() -> None:
    # 8 full repeats × 4 + 2 extra  →  CE 8.2
    assert ce_from_length(34, period=4, corr_value=0) == 8.2


def test_ce_from_length_th01_common_alleles() -> None:
    # TH01 alleles 6, 7, 8, 9, 10 must be exact integers with corr_value=0
    for allele in range(6, 11):
        assert ce_from_length(allele * 4, period=4, corr_value=0) == float(allele)


def test_ce_from_length_d22s1045_microvariant() -> None:
    # D22S1045 period=3: 16 full repeats × 3 + 1 extra  →  CE 16.1
    assert ce_from_length(49, period=3, corr_value=0) == 16.1


def test_ce_undefined_for_multi_motif() -> None:
    assert ce_from_length(60, period=-1, corr_value=0) is None


# ---------------------------------------------------------------------------
# ce_from_brackets (P4 — compound-motif CE)
# ---------------------------------------------------------------------------

def test_ce_from_brackets_d21s11_allele31() -> None:
    """D21S11 canonical allele 31: sum of all [TCTA]/[TCTG] bracket counts."""
    isfg = "[TCTA]6 [TCTG]5 [TCTA]3 ta [TCTA]3 tca [TCTA]2 tccata [TCTA]12"
    assert ce_from_brackets(isfg) == 31.0


def test_ce_from_brackets_d21s11_allele29() -> None:
    isfg = "[TCTA]4 [TCTG]6 [TCTA]3 ta [TCTA]3 tca [TCTA]2 tccata [TCTA]11"
    assert ce_from_brackets(isfg) == 29.0


def test_ce_from_brackets_single_copy_counts() -> None:
    """Bare uppercase motif tokens (single occurrences) count as 1."""
    # D8S1179-like: TCTA TCTG [TCTA]11
    isfg = "TCTA TCTG [TCTA]11"
    assert ce_from_brackets(isfg) == 13.0


def test_ce_from_brackets_lowercase_ignored() -> None:
    """Lowercase uncounted nucleotides (P3) must not be counted."""
    isfg = "agt [GGAA]10 gg"
    assert ce_from_brackets(isfg) == 10.0


def test_ce_from_brackets_pure_brackets() -> None:
    isfg = "[AGAT]9"
    assert ce_from_brackets(isfg) == 9.0


def test_ce_from_brackets_empty() -> None:
    assert ce_from_brackets("") is None


def test_ce_from_brackets_no_motif_found() -> None:
    """All-lowercase (no motif match) returns None."""
    assert ce_from_brackets("atcgatcgatcg") is None


# ---------------------------------------------------------------------------
# Strand-aware compress_isfg (P2)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Lowercase uncounted nucleotides (P3 — Phillips 2018 Rule 2)
# ---------------------------------------------------------------------------

def test_compress_isfg_lowercase_d13s317_like() -> None:
    """D13S317: trailing 'aat' must be lowercase."""
    # [TATC]11 + AAT (3 uncounted bases)
    seq = "TATC" * 11 + "AAT"
    out = compress_isfg(seq, motif="TATC")
    assert out == "[TATC]11 aat"


def test_compress_isfg_lowercase_d21s11_internal() -> None:
    """Internal uncounted tracts must be lowercase; bracket blocks stay uppercase."""
    # Simplified D21S11 fragment: [TCTA]3 TA [TCTA]2
    seq = "TCTA" * 3 + "TA" + "TCTA" * 2
    out = compress_isfg(seq, motif="TCTA,TCTG")
    assert "[TCTA]3" in out
    assert "[TCTA]2" in out
    assert " ta " in out  # internal uncounted must be lowercase


def test_compress_isfg_lowercase_d10s1248_flanks() -> None:
    """D10S1248: 'agt' prefix and 'gg' suffix must be lowercase."""
    # AGT [GGAA]4 GG
    seq = "AGT" + "GGAA" * 4 + "GG"
    out = compress_isfg(seq, motif="GGAA")
    assert out == "agt [GGAA]4 gg"


def test_compress_isfg_motif_stays_uppercase() -> None:
    """Single motif occurrences (not bracketed) must stay uppercase."""
    # TCTA (1×) + TCTG (1×) + [TCTA]3
    seq = "TCTA" + "TCTG" + "TCTA" * 3
    out = compress_isfg(seq, motif="TCTA,TCTG")
    # All tokens are motif matches — no lowercase expected
    assert out == "TCTA TCTG [TCTA]3"
    assert out == out.upper() or "[" in out  # no accidental lowercase on motif tokens


def test_compress_isfg_strand_plus_unchanged() -> None:
    """strand="+" (default) must behave identically to before."""
    assert compress_isfg("AGATAGATAGAT", motif="AGAT", strand="+") == "[AGAT]3"


def test_compress_isfg_strand_minus_d5s818_like() -> None:
    """D5S818: consensus on + strand is ATCT repeats; with strand="-" should give [AGAT]n."""
    # 9 × ATCT (= RC of AGAT)
    seq = "ATCT" * 9  # 36 bp, pure ATCT run
    out = compress_isfg(seq, motif="AGAT", strand="-")
    assert out == "[AGAT]9"


def test_compress_isfg_strand_minus_csf1po_like() -> None:
    """CSF1PO HG00113 allele 1: CTATCTATCT... with strand="-" and motif AGAT."""
    # From HG00113 CSV: CTATCTATCTATCTATCTATCTATCTATCTATCT (34 bp)
    seq = "CTATCTATCTATCTATCTATCTATCTATCTATCT"
    out = compress_isfg(seq, motif="AGAT", strand="-")
    # RC of seq = AGATAGATAGATAGATAGATAGATAGATAGATAG (34 bp)
    # = [AGAT]8 + AG (2 uncounted)
    assert "[AGAT]8" in out


def test_compress_isfg_strand_minus_with_flank() -> None:
    """Reverse-strand marker: uncounted flanking bases must appear in output."""
    # GAT + [AGAT]9 on RC strand = [AGAT]9 + ATC on + strand
    seq = "ATC" + "ATCT" * 9  # 3 flank + 36 repeat = 39 bp
    out = compress_isfg(seq, motif="AGAT", strand="-")
    assert "[AGAT]9" in out


def test_motif_repeat_summary_strand_minus() -> None:
    """motif_repeat_summary with strand="-" should count AGAT runs in D5S818-like seq."""
    seq = "ATCT" * 9
    s = motif_repeat_summary(seq, "AGAT", strand="-")
    assert "AGATx9" in s
