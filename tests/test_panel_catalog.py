"""Tests for the allele catalog model and JSON I/O (Phase 2.1)."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.errors import PanelError
from frontstr.panel.catalog import (
    AlleleCatalog,
    CatalogEntry,
    dump_catalog,
    load_catalog,
)
from frontstr.panel.seed_strseq import StrSeqRecord, build_catalog


def test_load_demo_catalog(demo_catalog_json: Path) -> None:
    cat = load_catalog(demo_catalog_json)
    assert cat.version == "demo-0.1"
    assert cat.reference_build == "GRCh38"
    assert set(cat.markers()) == {"D3S1358", "TH01"}
    assert len(cat.by_marker("D3S1358")) == 2
    assert cat.by_marker("nonexistent") == []


def test_d3s1358_isoalleles_distinct(demo_catalog_json: Path) -> None:
    """14a and 14b share CE 14 and length but differ in ISFG/suffix."""
    cat = load_catalog(demo_catalog_json)
    d3 = {e.suffix: e for e in cat.by_marker("D3S1358")}
    assert set(d3) == {"a", "b"}
    assert d3["a"].ce == 14 and d3["b"].ce == 14
    assert len(d3["a"].sequence) == len(d3["b"].sequence)  # same CE/length
    assert d3["a"].sequence != d3["b"].sequence  # different internal structure
    assert d3["a"].isfg != d3["b"].isfg


def test_load_missing_catalog(tmp_path: Path) -> None:
    with pytest.raises(PanelError, match="not found"):
        load_catalog(tmp_path / "nope.json")


def test_load_invalid_json(tmp_path: Path) -> None:
    p = tmp_path / "bad.json"
    p.write_text("{not valid json")
    with pytest.raises(PanelError, match="parse error"):
        load_catalog(p)


def test_entry_rejects_non_nucleotide_sequence() -> None:
    with pytest.raises(ValueError, match="ACGT"):
        CatalogEntry(marker="X", sequence="ACXT", isfg="[AC]1")


def test_entry_uppercases_sequence() -> None:
    e = CatalogEntry(marker="X", sequence="acgt", isfg="x")
    assert e.sequence == "ACGT"


def test_dump_roundtrip(tmp_path: Path) -> None:
    cat = AlleleCatalog(
        version="v1",
        entries=[
            CatalogEntry(marker="TH01", sequence="AATG", isfg="[AATG]1", ce=1, suffix="a"),
        ],
    )
    p = tmp_path / "out.json"
    dump_catalog(cat, p)
    back = load_catalog(p)
    assert back.version == "v1"
    assert back.by_marker("TH01")[0].sequence == "AATG"
    assert back.by_marker("TH01")[0].suffix == "a"


def test_build_catalog_dedups_marker_sequence() -> None:
    recs = [
        StrSeqRecord(marker="TH01", sequence="AATG", isfg="[AATG]1", ce=1, suffix="a"),
        StrSeqRecord(marker="TH01", sequence="aatg", isfg="[AATG]1", ce=1, suffix="a"),
        StrSeqRecord(marker="TH01", sequence="AATGAATG", isfg="[AATG]2", ce=2),
    ]
    cat = build_catalog(recs, version="strseq-test", description="t")
    th01 = cat.by_marker("TH01")
    assert len(th01) == 2  # case-folded duplicate collapsed
    assert cat.source == "STRSeq"
    assert {e.accession for e in th01} == {None}
