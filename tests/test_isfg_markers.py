"""Golden-file tests for per-marker ISFG bracket structures (P6).

Each test uses a reconstructed consensus derived from the HG00113 profile and
verifies that compress_isfg produces the canonical ISFG notation per Phillips
2018 Table 3 after P2 (RC), P3 (lowercase), and P6 (strand correction).
"""

from __future__ import annotations

from frontstr.interp.isfg import ce_from_brackets, compress_isfg

# ---------------------------------------------------------------------------
# D2S441 — strand "+" (TCTA motif on plus strand)
# ---------------------------------------------------------------------------


def test_d2s441_strand_is_plus_not_minus() -> None:
    """D2S441 motif is TCTA on the + strand — strand='+' must produce brackets."""
    # HG00113 allele 1: t [TCTA]10 ttta TCTA  (49 bp)
    consensus = "T" + "TCTA" * 10 + "TTTA" + "TCTA"
    assert len(consensus) == 49
    out = compress_isfg(consensus, motif="TCTA,TCAA", strand="+")
    assert "[TCTA]10" in out
    assert "ttta" in out  # internal intercalation must be lowercase
    assert out.startswith("t ")  # leading uncounted base must be lowercase


def test_d2s441_strand_minus_is_wrong() -> None:
    """Applying strand='-' to D2S441 must NOT produce TCTA brackets."""
    consensus = "T" + "TCTA" * 10 + "TTTA" + "TCTA"
    out = compress_isfg(consensus, motif="TCTA,TCAA", strand="-")
    assert "[TCTA]" not in out  # RC destroys the TCTA pattern


def test_d2s441_ce_from_brackets() -> None:
    """D2S441 HG00113 allele 1: 10+1=11 TCTA units → CE 11.0."""
    consensus = "T" + "TCTA" * 10 + "TTTA" + "TCTA"
    out = compress_isfg(consensus, motif="TCTA,TCAA", strand="+")
    assert ce_from_brackets(out) == 11.0


def test_d2s441_allele2_structure() -> None:
    """D2S441 allele 2: t [TCTA]8 tctg — TCTG is a single-base variant, not TCAA."""
    # HG00113 allele 2 original: T [TCTA]8 TCTG  (37 bp)
    consensus = "T" + "TCTA" * 8 + "TCTG"
    out = compress_isfg(consensus, motif="TCTA,TCAA", strand="+")
    assert "[TCTA]8" in out
    # TCTG is not in motif list → must appear lowercase
    assert "tctg" in out


# ---------------------------------------------------------------------------
# D1S1656 — strand "-" (TAGA motif on minus strand)
# ---------------------------------------------------------------------------


def test_d1s1656_strand_minus_produces_taga_brackets() -> None:
    """D1S1656 consensus on + strand is ATCT repeats; RC gives [TAGA]n."""
    # HG00113 allele 1: 11 ATCT on + strand + trailing A = 45 bp
    consensus = "ATCT" * 11 + "A"
    assert len(consensus) == 45
    out = compress_isfg(consensus, motif="TAGA,CAGA", strand="-")
    assert "[TAGA]11" in out
    assert out.endswith("t")  # trailing uncounted base must be lowercase


def test_d1s1656_strand_plus_fails() -> None:
    """Without RC, D1S1656 + strand shows no TAGA brackets."""
    consensus = "ATCT" * 11 + "A"
    out = compress_isfg(consensus, motif="TAGA,CAGA", strand="+")
    assert "[TAGA]" not in out


def test_d1s1656_ce_from_brackets() -> None:
    """D1S1656 HG00113 allele 1: 11 TAGA units → CE 11.0."""
    consensus = "ATCT" * 11 + "A"
    out = compress_isfg(consensus, motif="TAGA,CAGA", strand="-")
    assert ce_from_brackets(out) == 11.0


# ---------------------------------------------------------------------------
# D12S391 — strand "+" (AGAT,AGAC compound motif)
# ---------------------------------------------------------------------------


def test_d12s391_compound_motif() -> None:
    """D12S391 uses two motifs; both runs must be bracketed separately."""
    # HG00113: at [AGAT]6 [AGAC]8 + long non-motif trailer
    consensus = "AT" + "AGAT" * 6 + "AGAC" * 8 + "GAGAGGGGATTTATTAGA"
    out = compress_isfg(consensus, motif="AGAT,AGAC", strand="+")
    assert "[AGAT]6" in out
    assert "[AGAC]8" in out
    # Flanking bases must be lowercase
    assert out.startswith("at ")
    assert "gagaggggatttattaga" in out


# ---------------------------------------------------------------------------
# D18S51 — trailing 'ag' microvariant (P7)
# ---------------------------------------------------------------------------


def test_d18s51_trailing_ag_lowercase() -> None:
    """D18S51 canonical structure: [AGAA]n ag — trailing 2 bp must be lowercase."""
    # Allele 13: 13 × AGAA + AG = 54 bp (ref_end extended by 1 vs original YAML)
    consensus = "AGAA" * 13 + "AG"
    assert len(consensus) == 54
    out = compress_isfg(consensus, motif="AGAA")
    assert out == "[AGAA]13 ag"


def test_d18s51_ce_allele13() -> None:
    """With corr_value=2 and 54 bp → CE 13.0 (not 12.3)."""
    from frontstr.interp.isfg import ce_from_length

    assert ce_from_length(54, period=4, corr_value=2) == 13.0


def test_d18s51_ce_allele12() -> None:
    """Allele 12: 50 bp with corr_value=2 → CE 12.0."""
    from frontstr.interp.isfg import ce_from_length

    assert ce_from_length(50, period=4, corr_value=2) == 12.0


def test_d18s51_old_53bp_gives_truncated_output() -> None:
    """Regression: the old 53 bp extract (missing G) produced [AGAA]13 a (uppercase A)."""
    # After P3, non-motif is lowercase → was 'a' not 'ag'
    seq_old = "AGAA" * 13 + "A"
    out = compress_isfg(seq_old, motif="AGAA")
    assert out == "[AGAA]13 a"  # shows truncation — only fixed by coord change


def test_d12s391_ce_from_brackets() -> None:
    """D12S391 HG00113: 6+8=14 repeat units → CE 14.0."""
    consensus = "AT" + "AGAT" * 6 + "AGAC" * 8 + "GAGAGGGGATTTATTAGA"
    out = compress_isfg(consensus, motif="AGAT,AGAC", strand="+")
    assert ce_from_brackets(out) == 14.0
