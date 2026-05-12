"""Server-side SVG charts for the FRONTStr HTML report.

Charts are rendered as plain SVG strings — no JS, no fetches, no Vega-Lite
dependency. This keeps the report deterministic (byte-stable for the same
inputs) and printable to PDF without any JavaScript-aware engine.

The two charts that matter forensically:

- :func:`electropherogram_svg` — per-locus bar chart of integer read counts
  per cluster, colored by classification status, with the expected-stutter
  envelope drawn as a dashed line.
- :func:`coverage_bar_svg` — per-marker coverage overview for the QC page.

Style is intentionally MinKNOW-ish: thin axes, soft colors, no chart junk.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# Color palette (synced with payload status_chip codes)
STATUS_COLORS: dict[str, str] = {
    "allele": "#2e7d32",
    "inexact_allele": "#9e9d24",
    "artefact": "#f9a825",
    "stutter": "#90caf9",
    "noise": "#bdbdbd",
    "deletion": "#7b1fa2",
    "no_data": "#cfd8dc",
    "pending": "#e0e0e0",
}
EXPECTED_STUTTER_COLOR = "#455a64"
GRID_COLOR = "#e0e0e0"
AXIS_COLOR = "#90a4ae"
TEXT_COLOR = "#37474f"


@dataclass(slots=True)
class _ChartBox:
    width: int
    height: int
    margin_top: int = 20
    margin_right: int = 16
    margin_bottom: int = 36
    margin_left: int = 44

    @property
    def inner_w(self) -> int:
        return self.width - self.margin_left - self.margin_right

    @property
    def inner_h(self) -> int:
        return self.height - self.margin_top - self.margin_bottom


def electropherogram_svg(
    marker_result: dict[str, Any],
    *,
    width: int = 640,
    height: int = 260,
) -> str:
    """Render one locus's evidence clusters as a colored bar chart.

    Args:
        marker_result: One entry from ``payload["results"]`` (the JSON-shape,
            not a :class:`MarkerResult`).
        width: SVG outer width in CSS pixels.
        height: SVG outer height in CSS pixels.

    Returns:
        Complete ``<svg>...</svg>`` string. Returns an "empty locus" placeholder
        when there are no clusters.
    """
    alleles = list(marker_result.get("alleles", []))
    if not alleles:
        return _empty_chart_svg(width, height, "no reads at locus")

    box = _ChartBox(width=width, height=height)
    max_cov = max(int(a.get("n_reads_total", 0)) for a in alleles) or 1
    max_expected = max(float(a.get("expected_stutter", 0.0) or 0.0) for a in alleles)
    y_max = max(max_cov, int(max_expected) + 1)

    bar_slot = box.inner_w / len(alleles)
    bar_pad = max(2, int(bar_slot * 0.18))
    bar_w = max(6, int(bar_slot - 2 * bar_pad))

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" height="{height}" role="img" '
        f'aria-label="Electropherogram for {_xml(marker_result.get("marker_name", ""))}">'
    ]
    parts.append(_axes(box, y_max, y_label="reads"))

    for idx, a in enumerate(alleles):
        cov = int(a.get("n_reads_total", 0))
        status = str(a.get("status", "pending"))
        ce = a.get("ce")
        isfg = str(a.get("isfg") or "")
        consensus_len = int(a.get("length_bp", 0))
        expected = float(a.get("expected_stutter", 0.0) or 0.0)
        cx = box.margin_left + bar_slot * idx + bar_slot / 2
        bar_x = cx - bar_w / 2
        bar_h = (cov / y_max) * box.inner_h if y_max else 0
        bar_y = box.margin_top + box.inner_h - bar_h
        color = STATUS_COLORS.get(status, "#bdbdbd")
        label = f"CE{ce}" if ce is not None else f"{consensus_len}bp"
        tip = f"{label} | {status} | {cov} reads | ES={expected:.1f}"
        if isfg:
            tip += f" | {isfg}"
        parts.append(
            f'<rect x="{bar_x:.1f}" y="{bar_y:.1f}" width="{bar_w}" height="{bar_h:.1f}" '
            f'fill="{color}" rx="2" ry="2"><title>{_xml(tip)}</title></rect>'
        )
        if expected > 0:
            ex_h = min(expected, y_max) / y_max * box.inner_h
            ex_y = box.margin_top + box.inner_h - ex_h
            parts.append(
                f'<line x1="{bar_x - 2:.1f}" x2="{bar_x + bar_w + 2:.1f}" '
                f'y1="{ex_y:.1f}" y2="{ex_y:.1f}" '
                f'stroke="{EXPECTED_STUTTER_COLOR}" stroke-width="1.5" '
                f'stroke-dasharray="3 2"/>'
            )
        parts.append(
            f'<text x="{cx:.1f}" y="{box.margin_top + box.inner_h + 14:.1f}" '
            f'text-anchor="middle" font-size="11" fill="{TEXT_COLOR}">{_xml(label)}</text>'
        )
        parts.append(
            f'<text x="{cx:.1f}" y="{bar_y - 4:.1f}" text-anchor="middle" '
            f'font-size="11" fill="{TEXT_COLOR}" font-weight="600">{cov}</text>'
        )

    parts.append("</svg>")
    return "".join(parts)


def coverage_bar_svg(
    coverage_table: list[dict[str, Any]],
    *,
    width: int = 880,
    height_per_row: int = 22,
    margin_left: int = 110,
    floor: int = 30,
) -> str:
    """Horizontal bar chart of per-marker coverage with a dropout reference line."""
    if not coverage_table:
        return _empty_chart_svg(width, 120, "no markers")

    n = len(coverage_table)
    height = 40 + n * height_per_row
    margin_right = 20
    inner_w = width - margin_left - margin_right
    max_cov = max(int(r.get("coverage", 0)) for r in coverage_table) or 1

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" height="{height}" role="img" aria-label="Coverage per marker">'
    ]
    floor_x = margin_left + (min(floor, max_cov) / max_cov) * inner_w
    parts.append(
        f'<line x1="{floor_x:.1f}" x2="{floor_x:.1f}" y1="20" y2="{height - 16}" '
        f'stroke="#ef9a9a" stroke-width="1" stroke-dasharray="3 2"/>'
    )
    parts.append(
        f'<text x="{floor_x:.1f}" y="14" text-anchor="middle" '
        f'font-size="10" fill="#c62828">dropout floor ({floor}x)</text>'
    )

    for i, row in enumerate(coverage_table):
        marker = str(row.get("marker", ""))
        cov = int(row.get("coverage", 0))
        chip = str(row.get("chip", "ok"))
        y = 26 + i * height_per_row
        bar_w = (cov / max_cov) * inner_w if max_cov else 0
        color = {
            "ok": "#1976d2",
            "tri": "#7b1fa2",
            "mixture": "#c62828",
            "discordant": "#ef6c00",
            "no_data": "#9e9e9e",
        }.get(chip, "#1976d2")
        parts.append(
            f'<text x="{margin_left - 6:.1f}" y="{y + height_per_row / 2 + 3:.1f}" '
            f'text-anchor="end" font-size="11" fill="{TEXT_COLOR}">{_xml(marker)}</text>'
        )
        parts.append(
            f'<rect x="{margin_left}" y="{y}" width="{bar_w:.1f}" '
            f'height="{height_per_row - 6}" fill="{color}" rx="2" ry="2">'
            f'<title>{_xml(marker)}: {cov} reads ({chip})</title></rect>'
        )
        parts.append(
            f'<text x="{margin_left + bar_w + 4:.1f}" y="{y + height_per_row / 2 + 3:.1f}" '
            f'font-size="11" fill="{TEXT_COLOR}">{cov}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def haplotype_stack_svg(
    marker_result: dict[str, Any],
    *,
    width: int = 320,
    height: int = 200,
) -> str:
    """Stacked-bar of HP1 / HP2 / none reads for each called allele."""
    called = list(marker_result.get("alleles_called", []))
    if not called:
        return _empty_chart_svg(width, height, "nothing called")
    box = _ChartBox(width=width, height=height, margin_left=44, margin_right=12)
    y_max = max(int(a.get("n_reads_total", 0)) for a in called) or 1
    slot = box.inner_w / len(called)
    bar_w = max(14, int(slot * 0.7))
    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" height="{height}" role="img" aria-label="Haplotype split">'
    ]
    parts.append(_axes(box, y_max, y_label="reads"))
    for idx, a in enumerate(called):
        cx = box.margin_left + slot * idx + slot / 2
        bar_x = cx - bar_w / 2
        hp1 = int(a.get("n_reads_hp1", 0))
        hp2 = int(a.get("n_reads_hp2", 0))
        none = int(a.get("n_reads_hp_none", 0))
        y_floor = float(box.margin_top + box.inner_h)
        for value, color, label in (
            (none, "#cfd8dc", "no HP"),
            (hp2, "#9c27b0", "HP2"),
            (hp1, "#1976d2", "HP1"),
        ):
            if value <= 0:
                continue
            h = (value / y_max) * box.inner_h
            y_floor -= h
            parts.append(
                f'<rect x="{bar_x:.1f}" y="{y_floor:.1f}" width="{bar_w}" '
                f'height="{h:.1f}" fill="{color}" rx="1">'
                f'<title>{_xml(label)}: {value}</title></rect>'
            )
        ce = a.get("ce")
        label = f"CE{ce}" if ce is not None else f"#{idx}"
        parts.append(
            f'<text x="{cx:.1f}" y="{box.margin_top + box.inner_h + 14:.1f}" '
            f'text-anchor="middle" font-size="11" fill="{TEXT_COLOR}">{_xml(label)}</text>'
        )
    legend_y = box.margin_top + 4
    for i, (color, label) in enumerate(
        ((STATUS_COLORS["allele"], "HP1"), ("#9c27b0", "HP2"), ("#cfd8dc", "no HP")),
    ):
        cx = width - 70 + i * 24
        parts.append(
            f'<rect x="{cx}" y="{legend_y}" width="10" height="10" fill="{color}" rx="1"/>'
            f'<text x="{cx + 14}" y="{legend_y + 9}" font-size="10" '
            f'fill="{TEXT_COLOR}">{_xml(label)}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def _axes(box: _ChartBox, y_max: float, *, y_label: str) -> str:
    """Common axes + gridlines block."""
    parts: list[str] = []
    # Horizontal gridlines at 0%, 25%, 50%, 75%, 100% of y_max.
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = box.margin_top + box.inner_h * (1 - frac)
        value = round(y_max * frac)
        parts.append(
            f'<line x1="{box.margin_left}" x2="{box.margin_left + box.inner_w}" '
            f'y1="{y:.1f}" y2="{y:.1f}" stroke="{GRID_COLOR}" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{box.margin_left - 6}" y="{y + 3:.1f}" text-anchor="end" '
            f'font-size="10" fill="{TEXT_COLOR}">{value}</text>'
        )
    parts.append(
        f'<line x1="{box.margin_left}" x2="{box.margin_left}" '
        f'y1="{box.margin_top}" y2="{box.margin_top + box.inner_h}" '
        f'stroke="{AXIS_COLOR}" stroke-width="1"/>'
    )
    parts.append(
        f'<line x1="{box.margin_left}" x2="{box.margin_left + box.inner_w}" '
        f'y1="{box.margin_top + box.inner_h}" y2="{box.margin_top + box.inner_h}" '
        f'stroke="{AXIS_COLOR}" stroke-width="1"/>'
    )
    parts.append(
        f'<text x="6" y="{box.margin_top + box.inner_h / 2:.1f}" '
        f'transform="rotate(-90 6,{box.margin_top + box.inner_h / 2:.1f})" '
        f'font-size="11" fill="{TEXT_COLOR}">{_xml(y_label)}</text>'
    )
    return "".join(parts)


def _empty_chart_svg(width: int, height: int, message: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" height="{height}">'
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#fafafa" rx="4"/>'
        f'<text x="{width / 2:.1f}" y="{height / 2:.1f}" text-anchor="middle" '
        f'font-size="12" fill="#9e9e9e">{_xml(message)}</text>'
        f"</svg>"
    )


def _xml(value: str) -> str:
    """Minimal XML attribute / text escape."""
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
