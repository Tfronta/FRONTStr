"""Tests for the server-side SVG chart helpers."""

from __future__ import annotations

from xml.etree import ElementTree as ET

from frontstr.report.ngs_display import build_ngs_panel
from frontstr.report.svg_charts import (
    HP_SERIES,
    allele_coverage_svg,
    coverage_bar_svg,
    electropherogram_svg,
    haplotype_stack_svg,
)


def _allele(
    idx: int,
    ce: float,
    cov: int,
    status: str = "allele",
    expected: float = 0.0,
    consensus_suffix: str = "",
) -> dict:
    consensus = ("AATG" * int(ce)) + consensus_suffix
    return {
        "cluster_index": idx,
        "ce": ce,
        "n_reads_total": cov,
        "n_reads_hp1": cov // 2,
        "n_reads_hp2": cov - cov // 2,
        "n_reads_hp_none": 0,
        "length_bp": len(consensus),
        "consensus": consensus,
        "isfg": f"[AATG]{int(ce)}",
        "expected_stutter": expected,
        "status": status,
        "fraction": cov / 121,
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
        r
        for r in root.findall(".//{http://www.w3.org/2000/svg}rect")
        if r.get("class") == "ngs-segment"
    ]
    assert len(labeled) == 2


def test_allele_coverage_empty() -> None:
    svg = allele_coverage_svg({"marker_name": "X", "ngs_panel": {"chart_groups": []}})
    root = ET.fromstring(svg)
    text = root.find(".//{http://www.w3.org/2000/svg}text")
    assert text is not None
    assert "NGS panel" in (text.text or "")


def test_coverage_bar_includes_low_coverage_line() -> None:
    table = [
        {"marker": "M1", "coverage": 140, "called": 120, "chip": "ok"},
        {"marker": "M2", "coverage": 40, "called": 18, "chip": "ok"},
        {"marker": "M3", "coverage": 0, "called": 0, "chip": "no_data"},
    ]
    svg = coverage_bar_svg(table, floor=20)
    root = ET.fromstring(svg)
    lines = root.findall(".//{http://www.w3.org/2000/svg}line")
    dashed = [
        ln
        for ln in lines
        if "dasharray" in (ln.get("stroke-dasharray") or "") or ln.get("stroke-dasharray") == "3 2"
    ]
    assert dashed, "low coverage floor reference line missing"
    text_blob = "".join((t.text or "") for t in root.findall(".//{http://www.w3.org/2000/svg}text"))
    assert "low coverage floor" in text_blob
    assert "dropout" not in text_blob
    assert "M1" in text_blob and "M2" in text_blob


def test_coverage_bar_plots_supporting_reads_not_spanning() -> None:
    """The floor is measured on supporting reads, so the bar has to be too.

    M2 spans 40 reads but only 18 support its genotype, which is under the
    20-read floor. Plotting the spanning total would put its bar on the safe
    side of a line that the marker is in fact under.
    """
    table = [{"marker": "M2", "coverage": 40, "called": 18, "chip": "ok"}]
    root = ET.fromstring(coverage_bar_svg(table, floor=20))
    labels = [(t.text or "") for t in root.findall(".//{http://www.w3.org/2000/svg}text")]
    assert "18" in labels
    assert "40" not in labels


def test_coverage_bar_falls_back_for_payloads_without_called() -> None:
    """A run JSON serialized before ``called`` existed still draws its bars."""
    root = ET.fromstring(coverage_bar_svg([{"marker": "M1", "coverage": 33, "chip": "ok"}]))
    labels = [(t.text or "") for t in root.findall(".//{http://www.w3.org/2000/svg}text")]
    assert "33" in labels


def test_haplotype_stack_renders() -> None:
    marker = {
        "alleles_called": [_allele(0, 9.0, 30), _allele(1, 8.0, 25)],
    }
    svg = haplotype_stack_svg(marker)
    root = ET.fromstring(svg)
    rects = root.findall(".//{http://www.w3.org/2000/svg}rect")
    assert len(rects) >= 4


def _every_chart_svg() -> dict[str, str]:
    """One rendering of each chart, for invariants that hold across all of them."""
    marker = {
        "marker_name": "TPOX",
        "alleles_called": [_allele(0, 9.0, 30), _allele(1, 8.0, 25)],
    }
    marker["ngs_panel"] = build_ngs_panel(marker)
    return {
        "allele_coverage": allele_coverage_svg(marker),
        "haplotype_stack": haplotype_stack_svg(marker),
        "coverage_bar": coverage_bar_svg(
            [{"marker": "TPOX", "coverage": 40, "called": 18, "chip": "ok"}], floor=20
        ),
        "empty": allele_coverage_svg({"marker_name": "X"}),
    }


def test_chart_text_follows_the_stylesheet_not_a_literal() -> None:
    """No chart may hardcode a text colour.

    The report ships a light and a dark theme as custom properties, and these
    SVGs are inlined into it. A literal is therefore only ever right in one of
    the two: the axis and marker labels were ``#37474f``, a slate chosen for a
    white background, which on the dark theme rendered dark text on a dark
    panel and was unreadable. Every ``<text>`` fill has to defer to the
    stylesheet so both themes stay legible.
    """
    for name, svg in _every_chart_svg().items():
        texts = ET.fromstring(svg).findall(".//{http://www.w3.org/2000/svg}text")
        assert texts, f"{name} rendered no text at all"
        for t in texts:
            fill = t.get("fill") or ""
            assert fill.startswith("var(--"), (
                f"{name}: <text> fill {fill!r} is a literal; use a var(--…) token "
                f"so the dark theme stays readable"
            )


def test_chart_backgrounds_follow_the_stylesheet() -> None:
    """The placeholder panel must not paint a white box onto a dark page."""
    svg = allele_coverage_svg({"marker_name": "X"})
    rects = ET.fromstring(svg).findall(".//{http://www.w3.org/2000/svg}rect")
    assert rects, "empty-state chart drew no background"
    for r in rects:
        assert (r.get("fill") or "").startswith("var(--")


def _hp_marker() -> dict:
    """A phased heterozygote: one allele all HP1, the other all HP2."""
    a0 = _allele(0, 9.3, 10)
    a0["n_reads_hp1"] = 0
    a0["n_reads_hp2"] = 10
    a0["n_reads_hp_none"] = 0
    a1 = _allele(1, 7.0, 7)
    a1["n_reads_hp1"] = 7
    a1["n_reads_hp2"] = 0
    a1["n_reads_hp_none"] = 0
    return {"marker_name": "TH01", "alleles_called": [a0, a1]}


def _rects(svg: str) -> list:
    return ET.fromstring(svg).findall(".//{http://www.w3.org/2000/svg}rect")


def test_haplotype_legend_colours_match_the_bars() -> None:
    """A legend that names a colour the chart never draws is a false statement.

    HP1 was drawn blue in the bars and green in the legend, because the legend
    reused STATUS_COLORS["allele"]. On a phased heterozygote the chart contained
    no green at all while the key insisted HP1 was green.
    """
    svg = haplotype_stack_svg(_hp_marker())
    fills = {r.get("fill") for r in _rects(svg)}
    bar_fills = {colour for _key, colour, _label in HP_SERIES}
    # Every fill in the drawing is a declared haplotype colour, legend included.
    assert fills <= bar_fills, f"unexpected fill(s): {fills - bar_fills}"
    # The two haplotypes actually present are both drawn and both keyed.
    by_label = {label: colour for _k, colour, label in HP_SERIES}
    assert by_label["HP1"] in fills
    assert by_label["HP2"] in fills


def test_haplotype_legend_fits_and_does_not_overlap() -> None:
    """Swatch, label, swatch: each item has to clear the one before it."""
    width = 320
    svg = haplotype_stack_svg(_hp_marker(), width=width)
    root = ET.fromstring(svg)
    swatches = sorted(float(r.get("x") or 0) for r in _rects(svg) if r.get("y") == "2")
    assert len(swatches) == len(HP_SERIES), "one swatch per haplotype series"
    texts = sorted(
        (float(t.get("x") or 0), (t.text or ""))
        for t in root.findall(".//{http://www.w3.org/2000/svg}text")
        if t.get("y") == "11"
    )
    assert len(texts) == len(HP_SERIES)
    # Each label starts after its own swatch and ends before the next swatch.
    for i, (tx, label) in enumerate(texts):
        assert tx > swatches[i]
        end = tx + len(label) * 5.6
        if i + 1 < len(swatches):
            assert end <= swatches[i + 1], f"{label!r} runs under the next swatch"
        else:
            assert end <= width, f"{label!r} falls off the right edge"
