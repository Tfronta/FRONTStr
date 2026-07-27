"""NGS-style panel rows and stacked-bar grouping for HTML reports.

Rows mirror STRhub-style tables (allele · coverage · repeat · sequence) while
using FRONTStr evidence metrics: integer read counts per cluster consensus.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any

from frontstr.interp.stutter import find_motif_runs

_BRACKET_REPEAT = re.compile(r"\[([A-Za-z]+)\](\d+)")


def _xml_escape(s: str) -> str:
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def highlight_repeat_spans(consensus: str, motif_field: str) -> str:
    """Render sequence with STRhub-style flank (muted) + repeat (teal pill) spans."""
    motifs = [m for m in motif_field.split(",") if m]
    if not consensus or not motifs:
        return f'<span class="seq-flank">{_xml_escape(consensus)}</span>' if consensus else ""

    def _flank(s: str) -> str:
        return f'<span class="seq-flank">{_xml_escape(s)}</span>' if s else ""

    runs = find_motif_runs(consensus, motifs)
    if not runs:
        return _flank(consensus)

    runs_sorted = sorted(runs, key=lambda r: r.start)
    parts: list[str] = []
    pos = 0
    for run in runs_sorted:
        if run.start > pos:
            parts.append(_flank(consensus[pos : run.start]))
        chunk = consensus[run.start : run.end]
        parts.append(f'<span class="repeat-highlight">{_xml_escape(chunk)}</span>')
        pos = run.end
    if pos < len(consensus):
        parts.append(_flank(consensus[pos:]))
    return "".join(parts)


def _format_numeric_label(value: float, source: str) -> str:
    s = f"{value:.4f}".rstrip("0").rstrip(".")
    if source == "delta_only" and abs(value) > 1e-9:
        return f"Δ{s}"
    return s


def repeat_group_int(allele: dict[str, Any], period: int) -> int:
    """Repeat-count bin for chart X-axis and isoallele grouping."""
    an = allele.get("allele_numeric")
    if an is not None:
        return round(float(an))
    ce = allele.get("ce")
    if ce is not None:
        return round(float(ce))
    isfg = str(allele.get("isfg") or "")
    match = _BRACKET_REPEAT.search(isfg)
    if match:
        return int(match.group(2))
    length_bp = int(allele.get("length_bp") or 0)
    if period > 0:
        return length_bp // period
    return int(allele.get("cluster_index") or 0)


def _panel_source_alleles(marker: dict[str, Any]) -> list[dict[str, Any]]:
    called = list(marker.get("alleles_called") or [])
    if called:
        return called
    out: list[dict[str, Any]] = []
    for a in marker.get("alleles") or []:
        if str(a.get("status")) in ("allele", "inexact_allele"):
            out.append(a)
    return out


def build_ngs_panel(marker: dict[str, Any]) -> dict[str, Any]:
    """Derive table rows + chart_groups from one serialized marker dict."""
    marker_name = str(marker.get("marker_name") or "")
    motif = str(marker.get("system", {}).get("motif") or "")
    period = int(marker.get("system", {}).get("period") or 0)
    total_reads = int(marker.get("total_reads") or 0)

    subtitle = (
        "Allele haplotypes from ONT read clusters (FRONTStr evidence layer); "
        "coverage is read count per consensus sequence."
    )

    raw = _panel_source_alleles(marker)
    if not raw:
        return {
            "title": "Next-Generation Sequencing Analysis",
            "subtitle": subtitle,
            "rows": [],
            "chart_groups": [],
        }

    buckets: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for a in raw:
        rg = repeat_group_int(a, period)
        buckets[rg].append(a)

    rows: list[dict[str, Any]] = []
    chart_groups: list[dict[str, Any]] = []

    for rg in sorted(buckets.keys()):
        group_alleles = buckets[rg]
        by_seq: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for a in group_alleles:
            by_seq[str(a.get("consensus") or "")].append(a)

        distinct: list[dict[str, Any]] = []
        for lst in by_seq.values():
            distinct.append(max(lst, key=lambda x: int(x.get("n_reads_total") or 0)))

        canonical = max(
            distinct,
            key=lambda x: (
                int(x.get("n_reads_total") or 0),
                -int(x.get("cluster_index") or 0),
            ),
        )
        canonical_ci = int(canonical.get("cluster_index") or 0)

        group_rows_ordered = sorted(
            distinct,
            key=lambda x: (
                int(x.get("n_reads_total") or 0),
                int(x.get("cluster_index") or 0),
            ),
        )

        segments: list[dict[str, Any]] = []
        iso_counter = 0

        for idx_stack, a in enumerate(group_rows_ordered):
            cov = int(a.get("n_reads_total") or 0)
            consensus = str(a.get("consensus") or "")
            ci = int(a.get("cluster_index") or 0)
            row_id = f"{marker_name}:{ci}"
            is_canonical = ci == canonical_ci

            iso_ord = 0
            if not is_canonical:
                iso_counter += 1
                iso_ord = iso_counter

            frac = float(a.get("fraction") or 0.0)
            pct = round(frac * 100, 2) if total_reads > 0 else 0.0
            allele_num = a.get("allele_numeric")
            src = str(a.get("allele_numeric_source") or "")
            if allele_num is not None and src not in {"", "deletion", "unavailable"}:
                display = _format_numeric_label(float(allele_num), src)
            elif a.get("ce") is not None:
                display = _format_numeric_label(float(a["ce"]), "period_ce")
            else:
                display = str(rg)

            rows.append(
                {
                    "id": row_id,
                    "allele_display": display,
                    "repeat_group": rg,
                    "is_isoallele": not is_canonical,
                    "iso_ordinal": iso_ord,
                    "coverage_reads": cov,
                    "coverage_fraction_pct": pct,
                    "repeat_sequence": str(a.get("isfg") or ""),
                    "full_sequence": consensus,
                    "full_sequence_html": highlight_repeat_spans(consensus, motif),
                    "status": str(a.get("status") or ""),
                    "cluster_index": ci,
                    "stack_index": idx_stack,
                }
            )
            segments.append(
                {
                    "row_id": row_id,
                    "coverage_reads": cov,
                    "stack_index": idx_stack,
                    "is_isoallele": not is_canonical,
                }
            )

        chart_groups.append({"repeat_group": rg, "segments": segments})

    return {
        "title": "Next-Generation Sequencing Analysis",
        "subtitle": subtitle,
        "rows": rows,
        "chart_groups": chart_groups,
    }


def build_strhub_projection(payload: dict[str, Any]) -> dict[str, Any]:
    """Optional STRhub-shaped JSON slice for tooling / future import."""
    meta = payload.get("meta") or {}
    markers_out: list[dict[str, Any]] = []
    for r in payload.get("results") or []:
        panel = r.get("ngs_panel") or {}
        marker_rows: list[dict[str, Any]] = []
        for row in panel.get("rows") or []:
            marker_rows.append(
                {
                    "allele": row["allele_display"],
                    "allele_iso": row["is_isoallele"],
                    "iso_ordinal": row["iso_ordinal"],
                    "coverage_reads": row["coverage_reads"],
                    "coverage_fraction": row["coverage_fraction_pct"],
                    "repeat_sequence": row["repeat_sequence"],
                    "full_sequence": row["full_sequence"],
                }
            )
        markers_out.append({"locus": r["marker_name"], "rows": marker_rows})

    return {
        "schema": "strhub.ngs_panel/v1",
        "source": "frontstr",
        "app_version": meta.get("app_version", ""),
        "sample": meta.get("sample_name", ""),
        "markers": markers_out,
    }


__all__ = [
    "build_ngs_panel",
    "build_strhub_projection",
    "highlight_repeat_spans",
    "repeat_group_int",
]
