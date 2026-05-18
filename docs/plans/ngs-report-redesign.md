# FRONTStr HTML Report — NGS Panel Redesign Plan

**Version:** 1.0 · **Implemented:** yes (see repo)

This document defines the STRhub-inspired **layout** with FRONTStr-native **metrics**: per-allele **read coverage** (`n_reads_total`) and **full consensus sequence**, not CE peaks or LongTR PDP.

## Principles

1. No CE-style electropherogram metaphor — NGS allele coverage chart + 4-column table.
2. **Coverage (reads)** is integer read count per cluster; optional `% of locus` in tooltip only.
3. Isoalleles: adjacent table rows; **stacked bar** at one repeat-group tick on the chart.
4. Standalone HTML (inlined CSS/JS); chart/table sync via `data-row-id`.
5. **`strhub`** JSON projection is additive under the main payload.

## Implementation map

| Module | Role |
|--------|------|
| `frontstr/report/ngs_display.py` | `build_ngs_panel()`, repeat highlighting HTML |
| `frontstr/report/payload.py` | Embeds `ngs_panel` per marker + top-level `strhub` |
| `frontstr/report/svg_charts.py` | `allele_coverage_svg()` stacked teal bars |
| `frontstr/report/html.py` | Inlines `svg_allele_coverage` |
| `templates/run_report.html.j2` | NGS panel + forensic `<details>` |
| `static/app.js` | `initNgsPanels()` |
| `static/styles.css` | Panel layout, selection, repeat highlight |

## Resolved decisions

- Y-axis and table primary metric: **reads**.
- X-axis: repeat group integer (`int(round(ce))` or bracket parse).
- Table order within group: ascending reads (iso/smaller row above canonical when smaller).

See git history / tests for behavioral fixtures.
