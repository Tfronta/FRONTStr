"""Score a FRONTStr batch against published genotypes for the same samples.

**A development benchmark, not a step in calling a sample.** It only works on
the handful of public samples that have published genotypes; run FRONTStr on
your own data and there is nothing here to run. See :mod:`benchmark`.

Reads the ``batch_summary.csv`` that ``frontstr batch`` writes, joins it to the
three-technology sheet via :mod:`benchmark.truth`, and reports concordance per
marker plus a line-by-line discordance table.

Illumina is the only external anchor — longTR and STRspy are ONT callers on the
same reads, so agreeing with them is not independent confirmation. They are
still reported on every discordance, because that is what makes a disagreement
adjudicable: FRONTStr contradicting Illumina while matching both ONT callers is
a different finding from FRONTStr standing alone.

Usage::

    python -m benchmark.compare \\
        --summary out/batch_summary.csv \\
        --workbook ~/Desktop/1000GEN-ONT-Merged-Compar.xlsx \\
        --out out/concordance
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from benchmark.truth import (
    NO_EXTERNAL_TRUTH,
    TRUTH_SHEET,
    Genotype,
    TruthCall,
    canonical_genotype,
    canonical_sample_id,
    index_by_sample_marker,
    load_truth,
)

app = typer.Typer(add_completion=False, help=__doc__)
console = Console()

#: Columns the batch summary uses for bookkeeping; everything else is a marker.
NON_MARKER_COLUMNS = frozenset({"sample_id", "role", "status", "error"})


class Verdict(StrEnum):
    """Outcome of one (sample, marker) comparison against Illumina."""

    MATCH = "match"
    MISMATCH = "mismatch"
    NO_CALL = "no_call"  # FRONTStr declined or could not call the locus
    NO_TRUTH = "no_truth"  # Illumina has nothing to compare against


@dataclass(frozen=True, slots=True)
class Comparison:
    """One (sample, marker) row of the scored table."""

    sample_id: str
    marker: str
    verdict: Verdict
    frontstr: Genotype
    illumina: Genotype
    longtr: Genotype
    strspy: Genotype

    @property
    def agrees_with_longtr(self) -> bool:
        return bool(self.frontstr) and self.frontstr == self.longtr

    @property
    def agrees_with_strspy(self) -> bool:
        return bool(self.frontstr) and self.frontstr == self.strspy


def parse_summary(path: Path) -> dict[tuple[str, str], Genotype]:
    """Read ``batch_summary.csv`` into ``{(sample, marker): genotype}``.

    Only rows whose ``status`` is ``ok`` contribute; a sample that errored has
    empty marker cells anyway, and counting those as no-calls would blame the
    caller for what was a pipeline failure.
    """
    if not path.exists():
        raise FileNotFoundError(f"Batch summary not found: {path}")

    calls: dict[tuple[str, str], Genotype] = {}
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or "sample_id" not in reader.fieldnames:
            raise ValueError(f"{path} does not look like a batch summary (no sample_id column)")
        markers = [c for c in reader.fieldnames if c not in NON_MARKER_COLUMNS]

        for row in reader:
            if (row.get("status") or "").strip() != "ok":
                continue
            sample_id = canonical_sample_id(row["sample_id"])
            for marker in markers:
                cell = (row.get(marker) or "").strip()
                # "12.0,14.0" for a het, "13.0" for a hom, "?" where the CE
                # could not be derived, empty where nothing was called.
                alleles = [p for p in cell.split(",") if p]
                calls[(sample_id, marker)] = canonical_genotype(list(alleles))
    return calls


def compare(
    summary: dict[tuple[str, str], Genotype],
    truth_index: dict[tuple[str, str], dict[str, TruthCall]],
) -> list[Comparison]:
    """Join called genotypes to truth and assign a verdict to each locus.

    Only samples present on both sides are scored. Use :func:`unscored_samples`
    to report the ones that were called but have no truth row — silently
    dropping them makes a 5-sample batch look like a 2-sample validation.
    """
    called_samples = {sample for sample, _ in summary}
    rows: list[Comparison] = []
    for key in sorted(truth_index):
        sample_id, marker = key
        if sample_id not in called_samples:
            continue  # sample was not in this batch

        techs = truth_index[key]
        illumina = _genotype_of(techs, "ILLUMINA")
        longtr = _genotype_of(techs, "longTR")
        strspy = _genotype_of(techs, "STRspy")
        called = summary.get(key, ())

        if not called:
            verdict = Verdict.NO_CALL
        elif not illumina:
            verdict = Verdict.NO_TRUTH
        elif called == illumina:
            verdict = Verdict.MATCH
        else:
            verdict = Verdict.MISMATCH

        rows.append(
            Comparison(
                sample_id=sample_id,
                marker=marker,
                verdict=verdict,
                frontstr=called,
                illumina=illumina,
                longtr=longtr,
                strspy=strspy,
            )
        )
    return rows


def _genotype_of(techs: dict[str, TruthCall], technology: str) -> Genotype:
    """Genotype this technology reported, or ``()`` if it reported nothing."""
    call = techs.get(technology)
    return call.genotype if call is not None else ()


def unscored_samples(
    summary: dict[tuple[str, str], Genotype],
    truth_index: dict[tuple[str, str], dict[str, TruthCall]],
) -> list[str]:
    """Samples the batch called that the workbook has no row for."""
    with_truth = {sample for sample, _ in truth_index}
    return sorted({sample for sample, _ in summary} - with_truth)


def render_per_marker(rows: list[Comparison]) -> Table:
    """Build the per-marker concordance table."""
    table = Table(title="Concordance vs Illumina, per marker")
    table.add_column("Marker")
    table.add_column("Scored", justify="right")
    table.add_column("Match", justify="right")
    table.add_column("Mismatch", justify="right")
    table.add_column("No call", justify="right")
    table.add_column("No truth", justify="right")
    table.add_column("Concordance", justify="right")

    by_marker: dict[str, Counter[Verdict]] = {}
    for row in rows:
        by_marker.setdefault(row.marker, Counter())[row.verdict] += 1

    for marker in sorted(by_marker):
        counts = by_marker[marker]
        scored = counts[Verdict.MATCH] + counts[Verdict.MISMATCH]
        pct = f"{100 * counts[Verdict.MATCH] / scored:.1f}%" if scored else "—"
        note = "  (no Illumina)" if marker in NO_EXTERNAL_TRUTH else ""
        table.add_row(
            marker + note,
            str(scored),
            str(counts[Verdict.MATCH]),
            str(counts[Verdict.MISMATCH]),
            str(counts[Verdict.NO_CALL]),
            str(counts[Verdict.NO_TRUTH]),
            pct,
        )

    totals = Counter(row.verdict for row in rows)
    scored = totals[Verdict.MATCH] + totals[Verdict.MISMATCH]
    table.add_section()
    table.add_row(
        "[bold]TOTAL[/bold]",
        f"[bold]{scored}[/bold]",
        f"[bold]{totals[Verdict.MATCH]}[/bold]",
        f"[bold]{totals[Verdict.MISMATCH]}[/bold]",
        f"[bold]{totals[Verdict.NO_CALL]}[/bold]",
        f"[bold]{totals[Verdict.NO_TRUTH]}[/bold]",
        f"[bold]{100 * totals[Verdict.MATCH] / scored:.1f}%[/bold]" if scored else "—",
    )
    return table


def render_discordances(rows: list[Comparison], limit: int) -> Table:
    """Build the discordance table, with both ONT callers alongside."""
    table = Table(title=f"Discordances and no-calls (first {limit})")
    for column in (
        "Sample",
        "Marker",
        "Verdict",
        "FRONTStr",
        "Illumina",
        "longTR",
        "STRspy",
        "ONT",
    ):
        table.add_column(column)

    interesting = [r for r in rows if r.verdict in (Verdict.MISMATCH, Verdict.NO_CALL)]
    for row in interesting[:limit]:
        agree = [
            name
            for name, ok in (("longTR", row.agrees_with_longtr), ("STRspy", row.agrees_with_strspy))
            if ok
        ]
        table.add_row(
            row.sample_id,
            row.marker,
            row.verdict.value,
            _fmt(row.frontstr),
            _fmt(row.illumina),
            _fmt(row.longtr),
            _fmt(row.strspy),
            "+".join(agree) if agree else "—",
        )
    return table


def _fmt(genotype: Genotype) -> str:
    """Render a genotype for display; ``—`` when nothing was reported."""
    return "/".join(genotype) if genotype else "—"


def write_rows(rows: list[Comparison], path: Path) -> None:
    """Write every scored locus to CSV, for the audit trail."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(
            [
                "sample_id",
                "marker",
                "verdict",
                "frontstr",
                "illumina",
                "longtr",
                "strspy",
                "agrees_longtr",
                "agrees_strspy",
                "external_truth_available",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.sample_id,
                    row.marker,
                    row.verdict.value,
                    _fmt(row.frontstr),
                    _fmt(row.illumina),
                    _fmt(row.longtr),
                    _fmt(row.strspy),
                    int(row.agrees_with_longtr),
                    int(row.agrees_with_strspy),
                    int(row.marker not in NO_EXTERNAL_TRUTH),
                ]
            )


@app.command()
def main(
    summary: Annotated[Path, typer.Option("--summary", "-s", help="batch_summary.csv")],
    workbook: Annotated[Path, typer.Option("--workbook", "-w", help="Truth workbook .xlsx")],
    out: Annotated[Path | None, typer.Option("--out", "-o", help="Directory for the CSVs")] = None,
    sheet: Annotated[str, typer.Option("--sheet", help="Workbook sheet to score against")] = (
        TRUTH_SHEET
    ),
    limit: Annotated[int, typer.Option("--limit", help="Discordance rows to print")] = 40,
) -> None:
    """Score a FRONTStr batch against the external truth workbook."""
    calls = parse_summary(summary)
    converted: list[tuple[str, str, str, float, str]] = []
    truth = index_by_sample_marker(load_truth(workbook, sheet=sheet, converted=converted))
    rows = compare(calls, truth)

    if not rows:
        console.print("[red]No samples in the batch summary matched the workbook.[/red]")
        raise typer.Exit(code=1)

    samples = {row.sample_id for row in rows}
    console.print(f"Scored [bold]{len(samples)}[/bold] sample(s), {len(rows)} loci")
    skipped = unscored_samples(calls, truth)
    if skipped:
        console.print(
            f"[yellow]Not scored — no row in {sheet}:[/yellow] {', '.join(skipped)}"
            f"  ({len(skipped)} of {len(samples) + len(skipped)} called)"
        )
    if converted:
        # Named, not counted. These cells were scored as an allele other than
        # the one the sheet stores, and which ones matters more than how many.
        console.print(
            f"\n[yellow]{len(converted)} workbook cell(s) written in repeat units, "
            "read as ISFG:[/yellow]"
        )
        for sample_id, marker, tech, stored, isfg in converted:
            console.print(f"  {sample_id}  {marker}  {tech}  {stored:g} -> {isfg}")
    console.print()
    console.print(render_per_marker(rows))
    console.print()
    console.print(render_discordances(rows, limit))

    if out is not None:
        write_rows(rows, out / "concordance_loci.csv")
        console.print(f"\nWrote [bold]{out / 'concordance_loci.csv'}[/bold]")


if __name__ == "__main__":
    app()
