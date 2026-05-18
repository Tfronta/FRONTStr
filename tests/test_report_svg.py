"""Tests for the server-side SVG chart helpers."""

from __future__ import annotations

from xml.etree import ElementTree as ET

from frontstr.report.ngs_display import build_ngs_panel
from frontstr.report.svg_charts import (
    allele_coverage_svg,
    coverage_bar_svg,
    electropherogram_svg,
    haplotype_stack_svg,
)


def _allele(idx: int, ce: float, cov: int, status: str = "allele",
            expected: float = 0.0, consensus_suffix: str = "") -> dict:
    consensus = ("AATG" * int(ce)) + consensus_suffix
    return {
        "cluster_index": idx, "ce": ce, "n_reads_total": cov,
        "n_reads_hp1": cov // 2, "n_reads_hp2": cov - cov // 2,
        "n_reads_hp_none": 0, "length_bp": len(consensus),
        "consensus": consensus,
        "isfg": f"[AATG]{int(ce)}", "expected_stutter": expected,
        "status": status, "fraction": cov / 121,
    }


def _marker_payload(**extra):
    marker = {
        "marker_name": "TH01",
        "total_reads": 121,
        "system": {"motif": "AATG", "period": 4},
        "alleles": [
            _allele(0, 9.0, 60, "allele"),
            _allele(1, 8.0, 55, "allele"),
            _allele(2, 7.0, 5, "stutter", expected=5.5),
        ],
        "alleles_called": [
            _allele(0, 9.0, 60, "allele"),
            _allele(1, 8.0, 55, "allele"),
        ],
    }
    marker.update(extra)
    marker["ngs_panel"] = build_ngs_panel(marker)
    return marker


def test_allele_coverage_renders_valid_svg() -> None:
    marker = _marker_payload()
    svg = allele_coverage_svg(marker)
    root = ET.fromstring(svg)
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    rects = root.findall(".//{http://www.w3.org/2000/svg}rect")
    assert len(rects) >= 2
    assert "TH01" in (root.get("aria-label") or "")
    labeled = [r for r in rects if r.get("class") == "ngs-segment"]
    assert len(labeled) == 2
    # Bars ordered by repeat group ascending (CE 8 before CE 9).
    assert [r.get("data-row-id") for r in labeled] == ["TH01:1", "TH01:0"]


def test_electropherogram_deprecated_alias_matches() -> None:
    marker = _marker_payload()
    assert electropherogram_svg(marker) == allele_coverage_svg(marker)


def test_allele_coverage_stacked_isoallele_group() -> None:
    marker = {
        "marker_name": "vWA",
        "total_reads": 42,
        "system": {"motif": "TAGA", "period": 4},
        "alleles_called": [
            {
                "cluster_index": 0,
                "ce": 17.0,
                "consensus": "AAA" + "TAGA" * 17 + "CCCC",
                "isfg": "[TAGA]17",
                "n_reads_total": 19,
                "fraction": 19 / 42,
                "status": "allele",
                "length_bp": 75,
            },
            {
                "cluster_index": 1,
                "ce": 17.0,
                "consensus": "TTT" + "TAGA" * 17 + "CCCC",
                "isfg": "[TAGA]17",
                "n_reads_total": 23,
                "fraction": 23 / 42,
                "status": "allele",
                "length_bp": 75,
            },
        ],
        "alleles": [],
    }
    marker["ngs_panel"] = build_ngs_panel(marker)
    svg = allele_coverage_svg(marker)
    root = ET.fromstring(svg)
    labeled = [
        r for r in root.findall(".//{http://www.w3.org/2000/svg}rect")
        if r.get("class") == "ngs-segment"
    ]
    assert len(labeled) == 2


def test_allele_coverage_empty() -> None:
    svg = allele_coverage_svg({"marker_name": "X", "ngs_panel": {"chart_groups": []}})
    root = ET.fromstring(svg)
    text = root.find(".//{http://www.w3.org/2000/svg}text")
    assert text is not None
    assert "NGS panel" in (text.text or "")


def test_coverage_bar_includes_dropout_line() -> None:
    table = [
        {"marker": "M1", "coverage": 120, "chip": "ok"},
        {"marker": "M2", "coverage": 18, "chip": "ok"},
        {"marker": "M3", "coverage": 0, "chip": "no_data"},
    ]
    svg = coverage_bar_svg(table, floor=30)
    root = ET.fromstring(svg)
    lines = root.findall(".//{http://www.w3.org/2000/svg}line")
    dashed = [ln for ln in lines if "dasharray" in (ln.get("stroke-dasharray") or "") or
              ln.get("stroke-dasharray") == "3 2"]
    assert dashed, "dropout floor reference line missing"
    text_blob = "".join((t.text or "") for t in root.findall(".//{http://www.w3.org/2000/svg}text"))
    assert "dropout floor" in text_blob
    assert "M1" in text_blob and "M2" in text_blob


def test_haplotype_stack_renders() -> None:
    marker = {
        "alleles_called": [_allele(0, 9.0, 30), _allele(1, 8.0, 25)],
    }
    svg = haplotype_stack_svg(marker)
    root = ET.fromstring(svg)
    rects = root.findall(".//{http://www.w3.org/2000/svg}rect")
    assert len(rects) >= 4
