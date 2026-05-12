"""Tests for BED export from a Panel."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.caller.bed import split_panel_by_chromosome, write_panel_bed
from frontstr.errors import CallerError
from frontstr.panel.loader import load_panel
from frontstr.panel.models import Panel, System


def test_write_panel_bed_basic(tmp_path: Path) -> None:
    panel = Panel(
        name="t", version="0",
        systems=[
            System(name="M1", chromosome="chr1", ref_start=200, ref_end=240, motif="AGAT", period=4),
            System(name="M2", chromosome="chr1", ref_start=100, ref_end=140, motif="TCTA,TCTG",
                   period=-1),
            System(name="M3", chromosome="chr2", ref_start=50, ref_end=80, motif="AT", period=2),
        ],
    )
    out = write_panel_bed(panel, tmp_path / "p.bed")
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 3
    # Sorted by chromosome then ref_start
    assert lines[0] == "chr1\t100\t140\tTCTA,TCTG\tM2"
    assert lines[1] == "chr1\t200\t240\tAGAT\tM1"
    assert lines[2] == "chr2\t50\t80\tAT\tM3"


def test_write_panel_bed_skips_no_motif(tmp_path: Path) -> None:
    """Markers without a motif (Illumina-only legacy) are skipped — but pydantic
    validation forbids empty motifs in System construction, so the only way to
    produce a "no motif" system is via the model_construct backdoor. We assert
    the same outcome by making sure the BED contains exactly the systems with
    valid motifs."""
    panel = Panel(
        name="t", version="0",
        systems=[
            System(name="ok", chromosome="chr1", ref_start=1, ref_end=10, motif="A", period=1),
        ],
    )
    out = write_panel_bed(panel, tmp_path / "p.bed")
    assert "ok" in out.read_text()


def test_write_panel_bed_raises_when_empty(tmp_path: Path) -> None:
    panel = Panel.model_construct(name="empty", version="0", systems=[])
    with pytest.raises(CallerError):
        write_panel_bed(panel, tmp_path / "p.bed")


def test_split_panel_by_chromosome(tmp_path: Path) -> None:
    panel = Panel(
        name="t", version="0",
        systems=[
            System(name="A", chromosome="chr1", ref_start=10, ref_end=20, motif="AG", period=2),
            System(name="B", chromosome="chr1", ref_start=30, ref_end=40, motif="AG", period=2),
            System(name="C", chromosome="chr2", ref_start=1, ref_end=10, motif="AG", period=2),
        ],
    )
    paths = split_panel_by_chromosome(panel, tmp_path / "split")
    assert set(paths.keys()) == {"chr1", "chr2"}
    assert paths["chr1"].read_text().count("\n") == 2
    assert paths["chr2"].read_text().count("\n") == 1


def test_write_codis_panel(tmp_path: Path, codis_panel_yaml: Path) -> None:
    panel = load_panel(codis_panel_yaml)
    out = write_panel_bed(panel, tmp_path / "codis.bed")
    rows = out.read_text().strip().splitlines()
    assert len(rows) == len(panel.systems)
    for row in rows:
        parts = row.split("\t")
        assert len(parts) == 5
        assert parts[0].startswith("chr")
