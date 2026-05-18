"""HTML report layer.

Public surface:

- :func:`build_report` — write a self-contained ``.html`` report.
- :func:`serialize_run` — shared serializer used by HTML and exports.
- :class:`RunContext` — run-level metadata accepted by :func:`build_report`.
- :func:`allele_coverage_svg`, :func:`coverage_bar_svg`,
  :func:`haplotype_stack_svg` — server-side SVG chart helpers.
"""

from frontstr.report.html import build_report
from frontstr.report.payload import RunContext, serialize_run
from frontstr.report.svg_charts import (
    allele_coverage_svg,
    coverage_bar_svg,
    electropherogram_svg,
    haplotype_stack_svg,
)

__all__ = [
    "RunContext",
    "allele_coverage_svg",
    "build_report",
    "coverage_bar_svg",
    "electropherogram_svg",
    "haplotype_stack_svg",
    "serialize_run",
]
