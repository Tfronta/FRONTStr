"""Tests for panel YAML loading and model validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.errors import PanelError
from frontstr.panel.loader import load_panel
from frontstr.panel.models import Panel, System


def test_load_example_panel(codis_panel_yaml: Path) -> None:
    panel = load_panel(codis_panel_yaml)
    assert panel.name == "CODIS 20 + sex"
    assert panel.reference_build == "GRCh38"
    assert len(panel.systems) >= 20
    assert panel.by_name("TPOX") is not None
    tpox = panel.by_name("TPOX")
    assert tpox is not None
    assert tpox.allow_triallelic is True
    assert tpox.tri_balanced_thr == 0.5


def test_load_missing_panel(tmp_path: Path) -> None:
    with pytest.raises(PanelError, match="not found"):
        load_panel(tmp_path / "nope.yaml")


def test_load_invalid_yaml(tmp_path: Path) -> None:
    p = tmp_path / "bad.yaml"
    p.write_text(":\n  - not valid yaml: [\n")
    with pytest.raises(PanelError):
        load_panel(p)


def test_load_yaml_with_invalid_motif(tmp_path: Path) -> None:
    p = tmp_path / "bad.yaml"
    p.write_text(
        "name: x\nversion: '1'\nsystems:\n"
        "  - name: M\n    chromosome: chr1\n    ref_start: 1\n    ref_end: 10\n"
        "    motif: NNNN\n    period: 4\n"
    )
    with pytest.raises(PanelError):
        load_panel(p)


def test_panel_by_chromosome_grouping() -> None:
    panel = Panel(
        name="t", version="0",
        systems=[
            System(name="A", chromosome="chr1", ref_start=100, ref_end=120, motif="AG", period=2),
            System(name="B", chromosome="chr1", ref_start=50, ref_end=80, motif="AG", period=2),
            System(name="C", chromosome="chr2", ref_start=10, ref_end=40, motif="AG", period=2),
        ],
    )
    groups = panel.by_chromosome()
    assert list(groups.keys()) == ["chr1", "chr2"]
    assert [s.name for s in groups["chr1"]] == ["B", "A"]  # sorted by start


def test_panel_rejects_duplicate_names() -> None:
    with pytest.raises(Exception, match="duplicate"):
        Panel(
            name="t", version="0",
            systems=[
                System(name="X", chromosome="chr1", ref_start=1, ref_end=2, motif="AG", period=2),
                System(name="X", chromosome="chr1", ref_start=3, ref_end=4, motif="AG", period=2),
            ],
        )
