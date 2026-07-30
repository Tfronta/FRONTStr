"""Tests for catalog annotation of alleles (Phase 2.1 / plan §6.1)."""

from __future__ import annotations

from pathlib import Path

from frontstr.interp.catalog import annotate_alleles, annotate_with_catalog
from frontstr.interp.models import Allele, AlleleStatus
from frontstr.panel.catalog import AlleleCatalog, CatalogEntry, load_catalog
from frontstr.panel.models import System

D3 = System(
    name="D3S1358",
    chromosome="chr3",
    ref_start=1,
    ref_end=60,
    motif="TCTA,TCTG",
    period=4,
)

SEQ_14A = "TCTA" + "TCTG" * 1 + "TCTA" * 12
SEQ_14B = "TCTA" + "TCTG" * 2 + "TCTA" * 11


def _allele(seq: str, *, ce: float | None = 14.0) -> Allele:
    return Allele(
        cluster_index=0,
        consensus=seq,
        length_bp=len(seq),
        n_reads_total=20,
        n_reads_hp1=0,
        n_reads_hp2=0,
        n_reads_hp_none=20,
        n_forward=20,
        n_reverse=0,
        mean_qual=30.0,
        ce=ce,
        isfg="(live)",
        bp_diff=0,
        is_deletion=False,
        status=AlleleStatus.ALLELE,
    )


def _catalog() -> AlleleCatalog:
    return AlleleCatalog(
        version="t",
        entries=[
            CatalogEntry(
                marker="D3S1358",
                sequence=SEQ_14A,
                isfg="TCTA TCTG [TCTA]12",
                ce=14,
                suffix="a",
                source="STRSeq",
            ),
            CatalogEntry(
                marker="D3S1358",
                sequence=SEQ_14B,
                isfg="TCTA [TCTG]2 [TCTA]11",
                ce=14,
                suffix="b",
                source="STRSeq",
            ),
        ],
    )


def test_exact_match_adopts_isfg_and_suffix() -> None:
    a = _allele(SEQ_14B)
    m = annotate_with_catalog(a, D3, _catalog())
    assert m.matched and m.exact and m.distance == 0
    assert a.iso.suffix == "b"
    assert a.iso.match_type == "exact"
    assert a.iso.distance == 0
    assert a.iso.source == "STRSeq"
    assert a.isfg == "TCTA [TCTG]2 [TCTA]11"


def test_isoalleles_get_distinct_hits() -> None:
    """14a vs 14b resolve to distinct catalog hits."""
    cat = _catalog()
    a, b = _allele(SEQ_14A), _allele(SEQ_14B)
    annotate_with_catalog(a, D3, cat)
    annotate_with_catalog(b, D3, cat)
    assert a.iso.suffix == "a" and b.iso.suffix == "b"
    assert a.isfg != b.isfg


def test_approximate_match_records_approx_match_type() -> None:
    # One substitution mid-core (index 28, inside the TCTA run); gap-grouping keeps
    # the surrounding runs in one block so the edit shows up as distance 1.
    mutated = SEQ_14A[:28] + ("G" if SEQ_14A[28] != "G" else "C") + SEQ_14A[29:]
    a = _allele(mutated)
    m = annotate_with_catalog(a, D3, _catalog())
    assert m.matched and not m.exact and m.distance == 1
    assert a.iso.suffix == "a"
    assert a.iso.match_type == "approx"
    assert a.iso.distance == 1
    assert a.isfg == "TCTA TCTG [TCTA]12"


def test_core_extraction_strips_flanks_before_match() -> None:
    """A flank-wrapped consensus still matches the core-only catalog entry."""
    flanked = "gctgaaaagctccc" + SEQ_14B + "tgagacagggtcttgc"
    a = _allele(flanked.upper())
    m = annotate_with_catalog(a, D3, _catalog())
    assert m.exact and m.distance == 0
    assert a.iso.suffix == "b"


def test_no_hit_keeps_live_isfg() -> None:
    a = _allele("TCTA" * 20)  # length far outside the window
    m = annotate_with_catalog(a, D3, _catalog())
    assert not m.matched
    assert a.iso.suffix is None
    assert a.isfg == "(live)"


def test_distant_sequence_within_length_no_match() -> None:
    """Same length, but > eps edits → not adopted, live ISFG kept."""
    scrambled = "AAAA" + "CCCC" * 2 + "GGGG" * 11  # len == SEQ_14B, very different
    a = _allele(scrambled)
    m = annotate_with_catalog(a, D3, _catalog())
    assert not m.matched
    assert a.isfg == "(live)"


def test_deletion_allele_skipped() -> None:
    a = Allele(
        cluster_index=0,
        consensus="",
        length_bp=0,
        n_reads_total=5,
        n_reads_hp1=0,
        n_reads_hp2=0,
        n_reads_hp_none=5,
        n_forward=5,
        n_reverse=0,
        mean_qual=30.0,
        ce=None,
        isfg="",
        bp_diff=0,
        is_deletion=True,
        status=AlleleStatus.DELETION,
    )
    m = annotate_with_catalog(a, D3, _catalog())
    assert not m.matched


def test_annotate_alleles_noop_when_catalog_none() -> None:
    a = _allele(SEQ_14B)
    annotate_alleles([a], D3, None)
    assert a.iso.suffix is None
    assert a.isfg == "(live)"


def test_demo_catalog_th01_microvariant(demo_catalog_json: Path) -> None:
    cat = load_catalog(demo_catalog_json)
    th01 = System(name="TH01", chromosome="chr11", ref_start=1, ref_end=40, motif="AATG", period=4)
    seq_93 = "AATG" * 6 + "ATG" + "AATG" * 3
    a = _allele(seq_93, ce=9.3)
    m = annotate_with_catalog(a, th01, cat)
    assert m.exact
    assert a.isfg == "[AATG]6 ATG [AATG]3"
