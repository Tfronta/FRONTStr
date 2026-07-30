"""Cohort view: every sample at one marker, one block per marker.

The per-sample report answers "is this profile right". Across a hundred samples
that is the wrong question to have to ask a hundred times. The questions a
cohort raises are different:

- which samples failed, and at which loci;
- whether one locus fails across the whole cohort, which says something about
  the panel rather than about the samples;
- which samples are worth opening individually.

None of those can be answered by a hundred separate documents, so this view
pivots the profile: the marker is the block, the samples are the rows. Prior
art is the HipSTR web UI, which lays a multi-sample run out the same way. The
columns are FRONTStr's own, since the fields differ.

Every sample name links to that sample's individual report, which is where the
per-locus evidence lives. This view deliberately does not duplicate it.

Built from the run JSONs the batch already wrote, the same way the tidy export
is: the samples ran in worker processes and the files are on disk, so shipping
payloads back through pickle would be a cost paid for nothing.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from frontstr.errors import FrontstrError


@dataclass(slots=True)
class CohortSample:
    """One sample's identity and where its individual report lives."""

    sample_id: str
    #: Relative href to the per-sample report, or ``""`` when it was not
    #: written. A link that 404s is worse than no link.
    report_href: str = ""
    n_called: int = 0
    n_markers: int = 0
    n_flagged: int = 0

    @property
    def call_rate(self) -> float:
        return self.n_called / self.n_markers if self.n_markers else 0.0


@dataclass(slots=True)
class MarkerBlock:
    """One marker, with a row per sample.

    ``rows`` are the per-sample ``profile_rows`` entries verbatim, so this view
    and the single-sample report render the same numbers from the same source.
    """

    marker: str
    rows: list[dict[str, Any]] = field(default_factory=list)

    @property
    def n_samples(self) -> int:
        return len(self.rows)

    @property
    def n_called(self) -> int:
        return sum(1 for r in self.rows if r.get("allele1_ce_label"))

    @property
    def n_flagged(self) -> int:
        return sum(1 for r in self.rows if r.get("flags"))

    @property
    def call_rate(self) -> float:
        return self.n_called / self.n_samples if self.n_samples else 0.0

    @property
    def flag_counts(self) -> list[tuple[str, str, int]]:
        """``(short, code, n)`` per flag raised at this marker, commonest first.

        A marker whose every sample carries the same flag is the signal this
        view exists for: that is a statement about the panel or the assay, not
        about a hundred unlucky samples.
        """
        counts: Counter[tuple[str, str]] = Counter()
        for row in self.rows:
            for flag in row.get("flags", []):
                counts[(flag.get("short", "?"), flag.get("code", ""))] += 1
        return [(short, code, n) for (short, code), n in counts.most_common()]


def build_cohort_payload(
    payloads: list[dict[str, Any]],
    *,
    report_hrefs: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Pivot per-sample payloads into marker blocks.

    Args:
        payloads: Canonical run payloads, one per sample.
        report_hrefs: ``{sample_id: relative href}`` for the individual
            reports. Samples absent from it render without a link.

    Raises:
        FrontstrError: If no payload carries any profile rows.
    """
    hrefs = report_hrefs or {}

    samples: list[CohortSample] = []
    blocks: dict[str, MarkerBlock] = {}
    marker_order: list[str] = []

    for payload in payloads:
        sample_id = payload.get("meta", {}).get("sample_name", "")
        rows = payload.get("profile_rows", [])
        called = flagged = 0
        for row in rows:
            marker = row.get("marker", "")
            if not marker:
                continue
            if marker not in blocks:
                blocks[marker] = MarkerBlock(marker=marker)
                marker_order.append(marker)
            # The link rides on the row: the template renders rows, and looking
            # the sample up per cell would be a join in a loop.
            blocks[marker].rows.append({**row, "report_href": hrefs.get(sample_id, "")})
            if row.get("allele1_ce_label"):
                called += 1
            if row.get("flags"):
                flagged += 1
        samples.append(
            CohortSample(
                sample_id=sample_id,
                report_href=hrefs.get(sample_id, ""),
                n_called=called,
                n_markers=len(rows),
                n_flagged=flagged,
            )
        )

    if not marker_order:
        raise FrontstrError("No profile rows in any run JSON; nothing to build a cohort view from")

    # Panel order, not alphabetical: it is the order every other view uses and
    # the order a reader of forensic profiles already has in their head.
    ordered = [blocks[m] for m in marker_order]
    for block in ordered:
        block.rows.sort(key=lambda r: str(r.get("sample", "")))

    return {
        "samples": sorted(samples, key=lambda s: s.sample_id),
        "blocks": ordered,
        "summary": _cohort_summary(samples, ordered),
    }


def _cohort_summary(samples: list[CohortSample], blocks: list[MarkerBlock]) -> dict[str, Any]:
    """Headline numbers, plus the markers that are worst across the cohort."""
    total_cells = sum(b.n_samples for b in blocks)
    called_cells = sum(b.n_called for b in blocks)
    flagged_cells = sum(b.n_flagged for b in blocks)

    # Named rather than merely counted: "12% of calls are flagged" tells nobody
    # where to look, and the point of the cohort view is to point somewhere.
    worst_markers = sorted(blocks, key=lambda b: (b.call_rate, -b.n_flagged))[:5]
    worst_samples = sorted(samples, key=lambda s: (s.call_rate, -s.n_flagged))[:5]

    return {
        "n_samples": len(samples),
        "n_markers": len(blocks),
        "n_calls": total_cells,
        "n_called": called_cells,
        "n_flagged": flagged_cells,
        "call_rate": called_cells / total_cells if total_cells else 0.0,
        "worst_markers": [
            {"marker": b.marker, "call_rate": b.call_rate, "n_flagged": b.n_flagged}
            for b in worst_markers
            if b.call_rate < 1.0 or b.n_flagged
        ],
        "worst_samples": [
            {"sample": s.sample_id, "call_rate": s.call_rate, "n_flagged": s.n_flagged}
            for s in worst_samples
            if s.call_rate < 1.0 or s.n_flagged
        ],
    }


def build_cohort_report(
    payloads: list[dict[str, Any]],
    out_path: Path,
    *,
    report_hrefs: dict[str, str] | None = None,
    panel_name: str = "",
) -> Path:
    """Render the cohort view to ``out_path`` and return it."""
    from frontstr.report.html import render_template

    cohort = build_cohort_payload(payloads, report_hrefs=report_hrefs)
    rendered = render_template(
        "cohort_report.html.j2",
        {"cohort": cohort, "panel_name": panel_name},
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered, encoding="utf-8")
    return out_path
