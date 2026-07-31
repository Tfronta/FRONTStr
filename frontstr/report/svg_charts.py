"""Server-side SVG charts for the FRONTStr HTML report.

Charts are rendered as plain SVG strings — no JS, no fetches, no Vega-Lite
dependency. This keeps the report deterministic (byte-stable for the same
inputs) and printable to PDF without any JavaScript-aware engine.

Primary per-locus chart:

- :func:`allele_coverage_svg` — NGS-style stacked bars by repeat group with read
  counts per haplotype (isoalleles stack at one tick).

Also:

- :func:`coverage_bar_svg` — per-marker coverage overview for the QC page.

Style is intentionally MinKNOW-ish: thin axes, soft colors, no chart junk.
"""

from __future__ import annotations

import math
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
# Chrome colours follow the report's stylesheet rather than restating it. These
# charts are inlined into a document that defines both a light and a dark theme
# through custom properties, so a literal here is only ever right in one of
# them: the axis labels were #37474f, a slate picked for a white background,
# which on the dark theme was dark text on a dark panel and could not be read.
# The fallback after the comma is what a viewer without the stylesheet would
# use; inside the report the variable always wins.
GRID_COLOR = "var(--border, #e2e8f0)"
AXIS_COLOR = "var(--border-strong, #cbd5e1)"
TEXT_COLOR = "var(--text, #0f172a)"
MUTED_TEXT_COLOR = "var(--text-faint, #94a3b8)"
DANGER_COLOR = "var(--danger, #dc2626)"
SURFACE_COLOR = "var(--muted, #f1f5f9)"

TEAL_BAR = "#0d9488"
TEAL_BAR_ISO = "#14b8a6"

#: Haplotype fills, in the order they stack from the axis upward. One table, so
#: the bars and the legend cannot disagree: the legend used to draw HP1 with
#: ``STATUS_COLORS["allele"]``, an unrelated green, while the bars drew it blue,
#: and the chart contained no green at all.
#: Untagged reads are the muted one on purpose: they are the absence of a
#: haplotype, not a third one. #cfd8dc was near-white, which on the dark theme
#: made the smallest series the loudest thing in the chart.
HP_SERIES: tuple[tuple[str, str, str], ...] = (
    ("n_reads_hp_none", "#94a3b8", "no HP"),
    ("n_reads_hp2", "#9c27b0", "HP2"),
    ("n_reads_hp1", "#1976d2", "HP1"),
)


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


def allele_coverage_svg(
    marker_result: dict[str, Any],
    *,
    width: int = 640,
    height: int = 280,
) -> str:
    """NGS-style stacked bars: one tick per repeat group, Y = read counts."""
    panel = marker_result.get("ngs_panel") or {}
    groups = list(panel.get("chart_groups") or [])
    marker_name = str(marker_result.get("marker_name") or "")
    if not groups:
        return _empty_chart_svg(width, height, "no alleles in NGS panel")

    peak_stack = 1
    peak_seg = 1
    for grp in groups:
        segs = grp.get("segments") or []
        stack_sum = sum(int(s.get("coverage_reads") or 0) for s in segs)
        peak_stack = max(peak_stack, stack_sum)
        for s in segs:
            peak_seg = max(peak_seg, int(s.get("coverage_reads") or 0))
    y_max = float(_nice_y_axis_max(max(peak_stack, peak_seg)))

    box = _ChartBox(
        width=width,
        height=height,
        margin_top=max(14, min(22, height // 13)),
        margin_bottom=max(36, min(50, height // 5 + 8)),
        margin_left=44,
        margin_right=16,
    )
    n_groups = len(groups)
    slot_w = box.inner_w / max(n_groups, 1)

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" height="{height}" role="img" '
        f'aria-label="Allele coverage for {_xml(marker_name)}">'
    ]
    parts.append(
        _axes(
            box,
            y_max,
            y_label="Read coverage",
            grid_dashed=True,
            x_label="Allele",
        )
    )

    y_floor = float(box.margin_top + box.inner_h)

    for gi, grp in enumerate(groups):
        rg = grp["repeat_group"]
        segments = sorted(
            grp.get("segments") or [],
            key=lambda s: int(s.get("stack_index") or 0),
        )
        n_seg = len(segments)
        base_bw = slot_w * 0.52
        bar_w = min(
            slot_w * 0.9,
            base_bw + min(slot_w * 0.38, 14 * max(0, n_seg - 1)),
        )
        cx = box.margin_left + slot_w * gi + slot_w / 2
        bar_x = cx - bar_w / 2

        y_cursor = y_floor
        for si, seg in enumerate(segments):
            cov = int(seg.get("coverage_reads") or 0)
            h = (cov / y_max) * box.inner_h if y_max else 0.0
            y_cursor -= h
            iso = bool(seg.get("is_isoallele"))
            fill = TEAL_BAR_ISO if iso else TEAL_BAR
            opacity = "0.85" if iso else "1"
            row_id = str(seg.get("row_id") or "")
            reads_tip = f"{cov} reads"
            parts.append(
                f'<rect class="ngs-segment" tabindex="0" role="button" '
                f'data-row-id="{_xml(row_id)}" x="{bar_x:.1f}" y="{y_cursor:.1f}" '
                f'width="{bar_w:.1f}" height="{max(h, 0):.1f}" fill="{fill}" '
                f'opacity="{opacity}" rx="2" ry="2" stroke="#e2f8f5" stroke-width="0.6">'
                f"<title>{_xml(reads_tip)}</title></rect>"
            )
            if si < len(segments) - 1:
                parts.append(
                    f'<line x1="{bar_x:.1f}" x2="{bar_x + bar_w:.1f}" '
                    f'y1="{y_cursor:.1f}" y2="{y_cursor:.1f}" '
                    'stroke="#ffffff" stroke-width="1.2" opacity="0.85"/>'
                )

        tick = str(grp.get("label") or rg)
        parts.append(
            f'<text x="{cx:.1f}" y="{box.margin_top + box.inner_h + 26:.1f}" '
            f'text-anchor="middle" font-size="11" fill="{TEXT_COLOR}">{_xml(tick)}</text>'
        )

    parts.append("</svg>")
    return "".join(parts)


def electropherogram_svg(
    marker_result: dict[str, Any],
    **kw: Any,
) -> str:
    """Deprecated alias for :func:`allele_coverage_svg`."""
    return allele_coverage_svg(marker_result, **kw)


def _nice_y_axis_max(raw: float) -> float:
    """Round read-count axis up with modest headroom."""
    if raw <= 0:
        return 1.0
    headroom = raw * 1.06 + 1.0
    if headroom <= 12:
        return float(math.ceil(headroom))
    step = 10 ** math.floor(math.log10(headroom))
    nice_step = step if headroom / step <= 5 else step * 2
    return float(math.ceil(headroom / nice_step) * nice_step)


def _row_depth(row: dict[str, Any]) -> int:
    """Reads supporting the genotype for one coverage-table row.

    Falls back to ``coverage`` for payloads serialized before ``called`` was
    added, so an older run JSON still renders rather than drawing every bar at
    zero.
    """
    if "called" in row:
        return int(row.get("called") or 0)
    return int(row.get("coverage") or 0)


def coverage_bar_svg(
    coverage_table: list[dict[str, Any]],
    *,
    width: int = 880,
    height_per_row: int = 22,
    margin_left: int = 110,
    floor: int = 20,
) -> str:
    """Per-marker bar chart of the reads supporting the genotype.

    The bar is ``called``, not the window's spanning total, because ``floor`` is
    the ``LOW_COVERAGE`` threshold and that is measured against the supporting
    reads. Drawing the line across spanning depth compared two different
    quantities and put the line in the wrong place on every marker.
    """
    if not coverage_table:
        return _empty_chart_svg(width, 120, "no markers")

    n = len(coverage_table)
    height = 40 + n * height_per_row
    margin_right = 20
    inner_w = width - margin_left - margin_right
    max_cov = max(_row_depth(r) for r in coverage_table) or 1

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" height="{height}" role="img" aria-label="Coverage per marker">'
    ]
    floor_x = margin_left + (min(floor, max_cov) / max_cov) * inner_w
    parts.append(
        f'<line x1="{floor_x:.1f}" x2="{floor_x:.1f}" y1="20" y2="{height - 16}" '
        f'stroke="{DANGER_COLOR}" stroke-width="1" stroke-dasharray="3 2" opacity="0.7"/>'
    )
    parts.append(
        f'<text x="{floor_x:.1f}" y="14" text-anchor="middle" '
        f'font-size="10" fill="{DANGER_COLOR}">low coverage floor ({floor}x)</text>'
    )

    for i, row in enumerate(coverage_table):
        marker = str(row.get("marker", ""))
        cov = _row_depth(row)
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
            f"<title>{_xml(marker)}: {cov} read(s) supporting the genotype ({chip})</title></rect>"
        )
        parts.append(
            f'<text x="{margin_left + bar_w + 4:.1f}" y="{y + height_per_row / 2 + 3:.1f}" '
            f'font-size="11" fill="{TEXT_COLOR}">{cov}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


#: Legend geometry, in px at font-size 10.
_LEGEND_SWATCH = 8
_LEGEND_SWATCH_GAP = 5
_LEGEND_ITEM_GAP = 18
_LEGEND_CHAR_W = 5.6


def _hp_legend(width: int, margin_right: int) -> str:
    """Right-aligned haplotype key, laid out from the measured label widths.

    The previous version spaced items by a flat 24 px from ``width - 70``. At
    font-size 10 "HP1" is wider than that, so every label ran under the next
    swatch and "no HP" fell off the right edge. Sitting at ``margin_top + 4`` it
    was also inside the plot area, free to land on top of a bar; it now sits in
    the top margin, above the axis.
    """
    items = [(color, label) for _key, color, label in reversed(HP_SERIES)]
    widths = [
        _LEGEND_SWATCH + _LEGEND_SWATCH_GAP + len(label) * _LEGEND_CHAR_W for _c, label in items
    ]
    total = sum(widths) + _LEGEND_ITEM_GAP * (len(items) - 1)
    x = max(2.0, width - margin_right - total)

    parts: list[str] = []
    for (color, label), item_w in zip(items, widths, strict=True):
        parts.append(
            f'<rect x="{x:.1f}" y="3" width="{_LEGEND_SWATCH}" height="{_LEGEND_SWATCH}" '
            f'fill="{color}" rx="1"/>'
            f'<text x="{x + _LEGEND_SWATCH + _LEGEND_SWATCH_GAP:.1f}" y="10" font-size="10" '
            f'fill="{TEXT_COLOR}">{_xml(label)}</text>'
        )
        x += item_w + _LEGEND_ITEM_GAP
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
        y_floor = float(box.margin_top + box.inner_h)
        for key, color, label in HP_SERIES:
            value = int(a.get(key, 0))
            if value <= 0:
                continue
            h = (value / y_max) * box.inner_h
            y_floor -= h
            parts.append(
                f'<rect x="{bar_x:.1f}" y="{y_floor:.1f}" width="{bar_w}" '
                f'height="{h:.1f}" fill="{color}" rx="1">'
                f"<title>{_xml(label)}: {value}</title></rect>"
            )
        ce = a.get("ce")
        label = f"CE{ce}" if ce is not None else f"#{idx}"
        parts.append(
            f'<text x="{cx:.1f}" y="{box.margin_top + box.inner_h + 14:.1f}" '
            f'text-anchor="middle" font-size="11" fill="{TEXT_COLOR}">{_xml(label)}</text>'
        )
    parts.append(_hp_legend(width, box.margin_right))
    parts.append("</svg>")
    return "".join(parts)


def _axes(
    box: _ChartBox,
    y_max: float,
    *,
    y_label: str,
    grid_dashed: bool = False,
    x_label: str = "",
) -> str:
    """Common axes + gridlines block."""
    parts: list[str] = []
    dash = ' stroke-dasharray="2 4"' if grid_dashed else ""
    # Horizontal gridlines at 0%, 25%, 50%, 75%, 100% of y_max.
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = box.margin_top + box.inner_h * (1 - frac)
        value = round(y_max * frac)
        parts.append(
            f'<line x1="{box.margin_left}" x2="{box.margin_left + box.inner_w}" '
            f'y1="{y:.1f}" y2="{y:.1f}" stroke="{GRID_COLOR}" '
            f'stroke-width="1"{dash}/>'
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
    if x_label:
        cx = box.margin_left + box.inner_w / 2
        parts.append(
            f'<text x="{cx:.1f}" y="{box.height - 10:.1f}" text-anchor="middle" '
            f'font-size="11" fill="{TEXT_COLOR}">{_xml(x_label)}</text>'
        )
    return "".join(parts)


def _empty_chart_svg(width: int, height: int, message: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" height="{height}">'
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="{SURFACE_COLOR}" rx="4"/>'
        f'<text x="{width / 2:.1f}" y="{height / 2:.1f}" text-anchor="middle" '
        f'font-size="12" fill="{MUTED_TEXT_COLOR}">{_xml(message)}</text>'
        f"</svg>"
    )


def _xml(value: str) -> str:
    """Minimal XML attribute / text escape."""
    return (
        value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")
    )
