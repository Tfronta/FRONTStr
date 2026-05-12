"""Tests for the server-side SVG chart helpers."""

from __future__ import annotations

from xml.etree import ElementTree as ET

from frontstr.report.svg_charts import (
    coverage_bar_svg,
    electropherogram_svg,
    haplotype_stack_svg,
)


def _allele(idx: int, ce: float, cov: int, status: str = "allele",
            expected: float = 0.0) -> dict:
    return {
        "cluster_index": idx, "ce": ce, "n_reads_total": cov,
        "n_reads_hp1": cov // 2, "n_reads_hp2": cov - cov // 2,
        "n_reads_hp_none": 0, "length_bp": int(ce * 4),
        "isfg": f"[AGAT]{int(ce)}", "expected_stutter": expected,
        "status": status,
    }


def test_electropherogram_renders_valid_svg() -> None:
    marker = {
        "marker_name": "TH01",
        "alleles": [
            _allele(0, 9.0, 60, "allele"),
            _allele(1, 8.0, 55, "allele"),
            _allele(2, 7.0, 5, "stutter", expected=5.5),
        ],
    }
    svg = electropherogram_svg(marker)
    root = ET.fromstring(svg)
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    # 3 bars
    rects = root.findall(".//{http://www.w3.org/2000/svg}rect")
    assert len(rects) >= 3
    # The aria-label must contain the marker name
    assert "TH01" in (root.get("aria-label") or "")
    # Allele bars use the green color from STATUS_COLORS
    allele_rects = [r for r in rects if r.get("fill") == "#2e7d32"]
    assert len(allele_rects) == 2


def test_electropherogram_empty() -> None:
    svg = electropherogram_svg({"marker_name": "X", "alleles": []})
    root = ET.fromstring(svg)
    text = root.find(".//{http://www.w3.org/2000/svg}text")
    assert text is not None
    assert "no reads" in (text.text or "")


def test_coverage_bar_includes_dropout_line() -> None:
    table = [
        {"marker": "M1", "coverage": 120, "chip": "ok"},
        {"marker": "M2", "coverage": 18, "chip": "ok"},
        {"marker": "M3", "coverage": 0, "chip": "no_data"},
    ]
    svg = coverage_bar_svg(table, floor=30)
    root = ET.fromstring(svg)
    lines = root.findall(".//{http://www.w3.org/2000/svg}line")
    # Among them must be the dashed floor reference line
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
    # 2 alleles x up to 3 stack pieces + 3 legend swatches
    assert len(rects) >= 4
