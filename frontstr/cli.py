"""FRONTStr command-line interface.

Stages of the CLI are organised as Typer subcommands so each can be wired up
independently. Phase 1 ships the ``inspect`` subcommand (for debugging input
files) and a stub ``run`` that orchestrates the full pipeline.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Any

import typer
from rich.console import Console
from rich.table import Table

from frontstr.errors import FrontstrError
from frontstr.ingest import detect_input, validate_bam
from frontstr.interp.models import Allele
from frontstr.interp.qc import QcThresholds
from frontstr.version import __version__

app = typer.Typer(
    name="frontstr",
    help="Forensic Ranked Output for Nanopore Tandem STR profiling.",
    no_args_is_help=True,
    rich_markup_mode="rich",
)
console = Console()


def _interpret_allele_cell(a: Allele) -> str:
    """Format one genotype cell: canonical allele label … (reads).

    Reads :attr:`Allele.number_label` rather than re-deriving a number, so the
    CLI can never disagree with the report or the exports about what an allele
    is called.
    """
    return f"{a.number_label or '?'}({a.n_reads_total})"


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"[bold]FRONTStr[/bold] {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    _version: Annotated[
        bool | None,
        typer.Option(
            "--version",
            "-V",
            callback=_version_callback,
            is_eager=True,
            help="Show version and exit.",
        ),
    ] = None,
) -> None:
    """Forensic Ranked Output for Nanopore Tandem STR profiling."""


@app.command("inspect")
def inspect(
    path: Annotated[Path, typer.Argument(help="Input file to inspect.")],
    expected_build: Annotated[
        str | None,
        typer.Option("--build", help="Expected reference build for BAM/CRAM (e.g. GRCh38)."),
    ] = None,
) -> None:
    """Detect format of an input and validate it if it is a BAM/CRAM.

    This is a debugging helper. It does **not** run the pipeline; it just tells
    you what FRONTStr would do with the file.
    """
    try:
        info = detect_input(path)
    except FrontstrError as exc:
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=1) from exc

    t = Table(title=f"Input: {path.name}")
    t.add_column("Field", style="bold")
    t.add_column("Value")
    t.add_row("kind", info.kind)
    t.add_row("gzipped", "yes" if info.gzipped else "no")
    t.add_row("aligned", str(info.aligned) if info.aligned is not None else "n/a")
    t.add_row("size", f"{info.size_bytes / (1024 * 1024):.2f} MB")
    console.print(t)

    if info.kind in ("bam", "ubam", "cram") and info.aligned:
        try:
            report = validate_bam(path, expected_build=expected_build)
        except FrontstrError as exc:
            console.print(f"[red]validation error:[/red] {exc}")
            raise typer.Exit(code=2) from exc
        v = Table(title="BAM validation")
        v.add_column("Field", style="bold")
        v.add_column("Value")
        v.add_row("samples", ", ".join(report.sample_names) or "—")
        v.add_row("refs (n)", str(len(report.reference_names)))
        v.add_row("median MAPQ", str(report.median_mapq))
        v.add_row("records inspected", str(report.n_records_inspected))
        console.print(v)
        if report.warnings:
            console.print("\n[yellow]warnings:[/yellow]")
            for w in report.warnings:
                console.print(f"  • {w}")


@app.command("run")
def run(
    _input: Annotated[Path, typer.Option("--input", "-i", help="FASTQ/BAM/CRAM input.")],
    sample: Annotated[str, typer.Option("--sample", "-s", help="Sample identifier.")],
    panel: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    reference: Annotated[Path, typer.Option("--reference", "-r", help="Reference FASTA.")],
    out: Annotated[Path, typer.Option("--out", "-o", help="Output directory.")],
    platform: Annotated[
        str, typer.Option("--platform", help="Sequencing platform: ont|hifi.")
    ] = "ont",
) -> None:
    """Run the full FRONTStr pipeline on a single sample.

    NOTE: not yet implemented end-to-end. Phase 1 of the roadmap wires the
    layers (ingest → align/passthrough → LongTR → evidence → interp → report).
    """
    _ = _input, sample, panel, reference, out, platform
    console.print(
        "[yellow]not implemented yet[/yellow] — see ROADMAP.md (Phase 1).\n"
        "Use `frontstr inspect <path>` to validate inputs in the meantime."
    )
    raise typer.Exit(code=64)  # EX_USAGE


@app.command("batch")
def batch(
    manifest: Annotated[Path, typer.Option("--manifest", help="Batch manifest TSV.")],
    panel_path: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    out: Annotated[Path, typer.Option("--out", "-o", help="Output directory.")],
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", help="Reference FASTA (required for CRAM input)."),
    ] = None,
    formats: Annotated[
        str,
        typer.Option(
            "--formats",
            help=("Comma-separated list. One or more of: profile,evidence,seqs,json,html."),
        ),
    ] = "profile,evidence,seqs,json,html",
    workers: Annotated[
        int, typer.Option("--workers", "-j", help="Parallel processes (default: 1).")
    ] = 1,
    platform: Annotated[str, typer.Option("--platform")] = "ont",
    operator: Annotated[str | None, typer.Option("--operator")] = None,
    run_id: Annotated[str | None, typer.Option("--run-id")] = None,
    min_mapq: Annotated[int, typer.Option("--min-mapq")] = 20,
    identity: Annotated[float, typer.Option("--identity")] = 0.97,
    analytical_thresh: Annotated[float, typer.Option("--analytical-thresh")] = 0.02,
    calling_thresh: Annotated[float, typer.Option("--calling-thresh")] = 0.10,
) -> None:
    """Run FRONTStr on a multi-sample batch from a manifest TSV.

    The manifest is a tab-separated file with columns ``sample_id``, ``bam``
    and optionally ``role`` (sample|positive_ctrl|negative_ctrl|reagent_blank).
    Lines starting with ``#`` are ignored. Example::

        # FRONTStr batch manifest
        sample_id\\tbam\\trole
        HG00113\\t/data/HG00113.bam\\tsample
        CTRL_POS\\t/data/ctrl_pos.bam\\tpositive_ctrl

    Per-sample outputs are written under ``<out>/<sample_id>/``.
    A ``batch_summary.csv`` with all called CEs is written to ``<out>/``.
    """
    import os

    from frontstr.batch import parse_manifest, run_batch
    from frontstr.panel.loader import load_panel

    wanted = frozenset(f.strip().lower() for f in formats.split(",") if f.strip())
    unknown = wanted - {"profile", "evidence", "seqs", "json", "html"}
    if unknown:
        console.print(f"[red]Unknown formats:[/red] {sorted(unknown)}")
        raise typer.Exit(code=2)

    try:
        entries = parse_manifest(manifest)
        panel = load_panel(panel_path)
    except FrontstrError as exc:
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    n = len(entries)
    effective_workers = max(1, min(workers, n, os.cpu_count() or 1))
    console.print(
        f"[bold]FRONTStr batch[/bold] — {n} sample(s), "
        f"{effective_workers} worker(s), panel: [cyan]{panel.name}[/cyan]"
    )

    completed: list[str] = []

    def _tick(sample_id: str) -> None:
        completed.append(sample_id)
        console.print(f"  [{len(completed)}/{n}] {sample_id}")

    try:
        results = run_batch(
            entries=entries,
            panel=panel,
            out_dir=out,
            reference_fasta=reference,
            formats=wanted,
            min_mapq=min_mapq,
            identity=identity,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
            platform=platform,
            operator=operator,
            run_id=run_id,
            workers=effective_workers,
            progress_callback=_tick,
        )
    except FrontstrError as exc:
        console.print(f"[red]batch error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    ok = sum(1 for r in results if r.status == "ok")
    errors = [r for r in results if r.status == "error"]

    console.print(f"\n[green]✓ {ok}/{n} samples succeeded[/green]")
    if errors:
        console.print(f"[red]✗ {len(errors)} error(s):[/red]")
        for r in errors:
            first_line = r.error.splitlines()[0] if r.error else "unknown"
            console.print(f"  [red]{r.sample_id}[/red]: {first_line}")

    summary = out / "batch_summary.csv"
    if summary.exists():
        size_kb = summary.stat().st_size / 1024
        console.print(f"[green]wrote[/green] {summary}  ({size_kb:.1f} KB)")

    if errors:
        raise typer.Exit(code=1)


@app.command("call")
def call(
    bam: Annotated[Path, typer.Option("--bam", help="Indexed sample BAM.")],
    panel_path: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    reference: Annotated[
        Path, typer.Option("--reference", "-r", help="Reference FASTA (indexed).")
    ],
    out_dir: Annotated[Path, typer.Option("--out", "-o", help="Output directory.")],
    platform: Annotated[str, typer.Option("--platform", help="ont|hifi.")] = "ont",
    chrom: Annotated[
        str | None, typer.Option("--chrom", help="Restrict to a single chromosome.")
    ] = None,
    phased: Annotated[bool, typer.Option("--phased", help="Pass --phased-bam.")] = False,
    binary: Annotated[
        str | None, typer.Option("--longtr-bin", help="Path to LongTR binary.")
    ] = None,
    parse_only: Annotated[
        Path | None,
        typer.Option("--parse-only", help="Skip subprocess; just parse this VCF."),
    ] = None,
) -> None:
    """Run LongTR on a sample and pretty-print the parsed results.

    Use ``--parse-only PATH`` to skip the LongTR call entirely and just
    inspect a pre-existing VCF (handy for offline development).
    """
    from frontstr.caller import LongTRRunner, parse_longtr_vcf
    from frontstr.panel.loader import load_panel

    if parse_only is not None:
        try:
            results = parse_longtr_vcf(parse_only)
        except FrontstrError as exc:
            console.print(f"[red]parse error:[/red] {exc}")
            raise typer.Exit(code=2) from exc
    else:
        try:
            panel = load_panel(panel_path)
            runner = LongTRRunner(
                panel=panel,
                reference=reference,
                platform=platform,
                phased=phased,
                binary=binary,
            )
            run = runner.run(bam=bam, out_dir=out_dir, chrom=chrom)
            results = run.results
            console.print(f"[green]LongTR finished[/green] → {run.vcf}")
        except FrontstrError as exc:
            console.print(f"[red]caller error:[/red] {exc}")
            raise typer.Exit(code=2) from exc

    t = Table(title=f"LongTR results ({len(results)} loci)")
    t.add_column("Marker", style="bold")
    t.add_column("Position")
    t.add_column("Motif")
    t.add_column("GT")
    t.add_column("Q", justify="right")
    t.add_column("DP", justify="right")
    t.add_column("Alleles (bp_diff)")
    for r in results:
        sample = next(iter(r.samples.values()), None)
        gt = "/".join(str(i) for i in sample.gt_indices) if sample and sample.gt_indices else "."
        q = f"{sample.posterior:.2f}" if sample and sample.posterior is not None else "."
        dp = str(sample.depth) if sample else "."
        bp_diffs = ",".join(str(a.bp_diff) for a in r.alleles)
        t.add_row(r.marker_name, f"{r.chrom}:{r.pos}", r.motif or "?", gt, q, dp, bp_diffs)
    console.print(t)


@app.command("doctor")
def doctor(
    bam: Annotated[Path, typer.Option("--bam", help="Indexed sample BAM or CRAM.")],
    panel_path: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    min_mapq: Annotated[int, typer.Option("--min-mapq")] = 20,
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", help="Reference FASTA (required for CRAM input)."),
    ] = None,
) -> None:
    """Pre-flight sanity check: BAM ↔ panel compatibility.

    Runs every marker through the pileup and prints:
      - whether each chromosome exists in the BAM @SQ headers
      - read counts at each locus (with and without MAPQ filter)
      - a chromosome-naming hint if 'chr' prefix mismatches

    Run this BEFORE ``frontstr report`` to debug "no data" results.
    """
    import pysam

    from frontstr.evidence.consensus import poa_backend_name
    from frontstr.evidence.pileup import pileup_locus
    from frontstr.panel.loader import load_panel

    backend = poa_backend_name()
    if backend:
        console.print(f"[green]POA backend:[/green] {backend}")
    else:
        console.print(
            "[red]POA backend: none[/red] — cluster consensus will fall back to a "
            "single unpolished read, degrading ISFG strings and iso-allele calls.\n"
            "  Fix: [bold]pip install 'frontstr[poa]'[/bold]"
        )

    try:
        panel = load_panel(panel_path)
    except FrontstrError as exc:
        console.print(f"[red]panel error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    is_cram = bam.suffix.lower() == ".cram"
    if is_cram and reference is None:
        console.print("[red]error:[/red] CRAM input requires --reference <fasta>")
        raise typer.Exit(code=2)
    open_kwargs: dict[str, Any] = {"reference_filename": str(reference)} if is_cram else {}
    try:
        af = pysam.AlignmentFile(str(bam), "rc" if is_cram else "rb", **open_kwargs)
    except (OSError, ValueError) as exc:
        console.print(f"[red]BAM error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    bam_chroms = set(af.references)
    panel_chroms = {s.chromosome for s in panel.systems}
    missing = panel_chroms - bam_chroms
    if missing:
        with_chr = {c if c.startswith("chr") else f"chr{c}" for c in bam_chroms}
        without_chr = {c.removeprefix("chr") for c in bam_chroms}
        if missing.issubset(with_chr):
            console.print(
                "[yellow]Hint:[/yellow] BAM has no 'chr' prefix; panel does. "
                "Use a panel without 'chr' or rename BAM @SQ headers."
            )
        elif missing.issubset({f"chr{c}" for c in without_chr}):
            console.print("[yellow]Hint:[/yellow] BAM has 'chr' prefix; panel does not.")

    t = Table(title=f"FRONTStr doctor — {panel.name} vs {bam.name}")
    t.add_column("Marker", style="bold")
    t.add_column("Chrom")
    t.add_column("Region")
    t.add_column("@SQ?")
    t.add_column("Reads @MAPQ0", justify="right")
    t.add_column(f"Reads @MAPQ{min_mapq}", justify="right")
    t.add_column("Status")
    for s in panel.systems:
        chrom_ok = s.chromosome in bam_chroms
        n_raw = 0
        n_filtered = 0
        if chrom_ok:
            try:
                obs_raw = pileup_locus(
                    bam,
                    s.chromosome,
                    s.ref_start - 1,
                    s.ref_end,
                    min_mapq=0,
                    reference_fasta=reference,
                )
                obs_filt = pileup_locus(
                    bam,
                    s.chromosome,
                    s.ref_start - 1,
                    s.ref_end,
                    min_mapq=min_mapq,
                    reference_fasta=reference,
                )
                n_raw, n_filtered = len(obs_raw), len(obs_filt)
            except FrontstrError as exc:
                t.add_row(
                    s.name,
                    s.chromosome,
                    f"{s.ref_start}-{s.ref_end}",
                    "yes",
                    "?",
                    "?",
                    f"[red]{exc}[/red]",
                )
                continue
        status: str
        if not chrom_ok:
            status = "[red]chrom missing[/red]"
        elif n_filtered == 0 and n_raw > 0:
            status = "[yellow]all reads below MAPQ[/yellow]"
        elif n_filtered == 0:
            status = "[red]no reads[/red]"
        elif n_filtered < 10:
            status = "[yellow]low cov[/yellow]"
        else:
            status = "[green]ok[/green]"
        t.add_row(
            s.name,
            s.chromosome,
            f"{s.ref_start}-{s.ref_end}",
            "yes" if chrom_ok else "no",
            str(n_raw),
            str(n_filtered),
            status,
        )
    console.print(t)
    af.close()


@app.command("export")
def export_cmd(
    bam: Annotated[Path, typer.Option("--bam", help="Indexed sample BAM or CRAM.")],
    panel_path: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    out_dir: Annotated[Path, typer.Option("--out-dir", "-o", help="Output directory.")],
    formats: Annotated[
        str,
        typer.Option(
            "--formats",
            help=(
                "Comma-separated list. One or more of: "
                "profile,evidence,seqs,json,json-compact,html,vcf,xlsx. "
                "vcf requires --reference."
            ),
        ),
    ] = "profile,evidence,seqs,json",
    sample_name: Annotated[
        str | None, typer.Option("--sample", help="Sample name (defaults to BAM stem).")
    ] = None,
    longtr_vcf: Annotated[
        Path | None, typer.Option("--longtr-vcf", help="Optional LongTR VCF.")
    ] = None,
    operator: Annotated[str | None, typer.Option("--operator")] = None,
    run_id: Annotated[str | None, typer.Option("--run-id")] = None,
    platform: Annotated[str, typer.Option("--platform")] = "ont",
    min_mapq: Annotated[int, typer.Option("--min-mapq")] = 20,
    identity: Annotated[float, typer.Option("--identity")] = 0.97,
    analytical_thresh: Annotated[float, typer.Option("--analytical-thresh")] = 0.02,
    calling_thresh: Annotated[float, typer.Option("--calling-thresh")] = 0.10,
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", help="Reference FASTA (required for CRAM input)."),
    ] = None,
    low_coverage_reads: Annotated[
        int,
        typer.Option(
            "--low-coverage-reads",
            help="Called loci below this many reads raise LOW_COVERAGE.",
        ),
    ] = QcThresholds.model_fields["low_coverage_reads"].default,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Log every marker call, not just the run.")
    ] = False,
) -> None:
    """Run the full pipeline and write one or more export files.

    A JSONL process log is always written to ``<out-dir>/frontstr.log.jsonl``
    alongside the exports; the canonical JSON carries the audit record.

    Example::

        frontstr export --bam s.bam --panel p.yaml -o exports/ \
            --formats profile,evidence,seqs,json,html
    """
    import logging

    from frontstr.caller import parse_longtr_vcf
    from frontstr.exports import (
        write_evidence_csv,
        write_profile_csv,
        write_run_json,
        write_run_vcf,
        write_run_xlsx,
        write_seqs_csv,
    )
    from frontstr.interp import index_longtr_results, interpret_run
    from frontstr.log import PROCESS_LOG_NAME, configure_logging
    from frontstr.panel.loader import load_panel
    from frontstr.report import RunContext, build_report, serialize_run

    wanted = {f.strip().lower() for f in formats.split(",") if f.strip()}
    unknown = wanted - {
        "profile",
        "evidence",
        "seqs",
        "json",
        "json-compact",
        "html",
        "vcf",
        "xlsx",
    }
    if unknown:
        console.print(f"[red]Unknown formats:[/red] {sorted(unknown)}")
        raise typer.Exit(code=2)
    if "vcf" in wanted and reference is None:
        # Fail before the pipeline runs rather than after several minutes of work.
        console.print(
            "[red]error:[/red] --formats vcf needs --reference <fasta>. REF must be "
            "the real reference sequence; there is no meaningful placeholder."
        )
        raise typer.Exit(code=2)

    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / PROCESS_LOG_NAME
    configure_logging(log_path, level=logging.DEBUG if verbose else logging.INFO)
    qc_thresholds = QcThresholds(low_coverage_reads=low_coverage_reads)

    try:
        panel = load_panel(panel_path)
        longtr_map = index_longtr_results(parse_longtr_vcf(longtr_vcf)) if longtr_vcf else None
        results = interpret_run(
            bam=bam,
            panel=panel,
            longtr_results=longtr_map,
            min_mapq=min_mapq,
            identity_threshold=identity,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
            reference_fasta=reference,
            qc_thresholds=qc_thresholds,
        )
        context = RunContext(
            sample_name=sample_name or bam.stem,
            panel_name=panel.name,
            panel_version=panel.version,
            bam_path=bam,
            longtr_vcf_path=longtr_vcf,
            platform=platform,
            operator=operator,
            run_id=run_id,
            reference_build=panel.reference_build,
            qc_thresholds=qc_thresholds,
        )
        payload = serialize_run(results, context)
    except FrontstrError as exc:
        console.print(f"[red]export error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    stem = context.sample_name
    written: list[Path] = []
    if "profile" in wanted:
        written.append(write_profile_csv(payload, out_dir / f"{stem}.profile.csv"))
    if "evidence" in wanted:
        written.append(write_evidence_csv(payload, out_dir / f"{stem}.evidence.csv"))
    if "seqs" in wanted:
        written.append(write_seqs_csv(payload, out_dir / f"{stem}.seqs.csv"))
    if "json" in wanted:
        written.append(write_run_json(payload, out_dir / f"{stem}.json", mode="pretty"))
    if "json-compact" in wanted:
        written.append(write_run_json(payload, out_dir / f"{stem}.min.json", mode="compact"))
    if "vcf" in wanted:
        written.append(write_run_vcf(payload, out_dir / f"{stem}.vcf", reference_fasta=reference))
    if "xlsx" in wanted:
        written.append(write_run_xlsx(payload, out_dir / f"{stem}.xlsx"))
    if "html" in wanted:
        written.append(build_report(results, context, out_dir / f"{stem}.html"))

    for path in written:
        size_kb = path.stat().st_size / 1024
        console.print(f"[green]wrote[/green] {path}  ({size_kb:.1f} KB)")

    audit = payload["audit"]
    sev = audit["severity_counts"]
    console.print(f"[green]wrote[/green] {log_path}  (process log)")
    if sev.get("error") or sev.get("warn"):
        console.print(
            f"[yellow]review:[/yellow] {sev.get('error', 0)} error / "
            f"{sev.get('warn', 0)} warning flag(s) across "
            f"{len(audit['markers_needing_review'])} marker(s): "
            f"{', '.join(audit['markers_needing_review'])}"
        )
    else:
        console.print("[green]no review flags raised[/green]")


@app.command("report")
def report(
    bam: Annotated[Path, typer.Option("--bam", help="Indexed sample BAM or CRAM.")],
    panel_path: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    out: Annotated[Path, typer.Option("--out", "-o", help="Output HTML path.")],
    sample_name: Annotated[
        str | None, typer.Option("--sample", help="Sample name (defaults to BAM stem).")
    ] = None,
    longtr_vcf: Annotated[
        Path | None, typer.Option("--longtr-vcf", help="Optional LongTR VCF.")
    ] = None,
    operator: Annotated[str | None, typer.Option("--operator")] = None,
    run_id: Annotated[str | None, typer.Option("--run-id")] = None,
    platform: Annotated[str, typer.Option("--platform")] = "ont",
    min_mapq: Annotated[int, typer.Option("--min-mapq")] = 20,
    identity: Annotated[float, typer.Option("--identity")] = 0.97,
    analytical_thresh: Annotated[float, typer.Option("--analytical-thresh")] = 0.02,
    calling_thresh: Annotated[float, typer.Option("--calling-thresh")] = 0.10,
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", help="Reference FASTA (required for CRAM input)."),
    ] = None,
) -> None:
    """Run the full pipeline and emit a self-contained HTML report.

    Usage::

        frontstr report --bam s.bam --panel codis.yaml --out s.html [--longtr-vcf s.vcf]
    """
    from frontstr.caller import parse_longtr_vcf
    from frontstr.interp import index_longtr_results, interpret_run
    from frontstr.panel.loader import load_panel
    from frontstr.report import RunContext, build_report

    try:
        panel = load_panel(panel_path)
        longtr_map = index_longtr_results(parse_longtr_vcf(longtr_vcf)) if longtr_vcf else None
        results = interpret_run(
            bam=bam,
            panel=panel,
            longtr_results=longtr_map,
            min_mapq=min_mapq,
            identity_threshold=identity,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
            reference_fasta=reference,
        )
        context = RunContext(
            sample_name=sample_name or bam.stem,
            panel_name=panel.name,
            panel_version=panel.version,
            bam_path=bam,
            longtr_vcf_path=longtr_vcf,
            platform=platform,
            operator=operator,
            run_id=run_id,
            reference_build=panel.reference_build,
        )
        out_path = build_report(results, context, out)
    except FrontstrError as exc:
        console.print(f"[red]report error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    size_kb = out_path.stat().st_size / 1024
    console.print(f"[green]✓ Report written:[/green] {out_path} ({size_kb:.1f} KB)")


@app.command("tidy")
def tidy(
    out_dir: Annotated[Path, typer.Option("--out-dir", "-o", help="Where to write the dataset.")],
    json_paths: Annotated[
        list[Path] | None,
        typer.Option("--json", help="Run JSON. Repeat, or use --from-dir."),
    ] = None,
    from_dir: Annotated[
        Path | None,
        typer.Option("--from-dir", help="Directory to search recursively for run JSONs."),
    ] = None,
    stem: Annotated[str, typer.Option("--stem", help="Output filename stem.")] = "cohort_tidy",
) -> None:
    """Build a cohort-scale tidy dataset from run JSONs.

    One row per sample x marker x allele — the shape a concordance study wants.
    Written as CSV and, when pyarrow is available, Parquet.

    Built from the canonical JSONs rather than by re-running, so a dataset can
    be rebuilt at any time, and runs from different batches combined::

        frontstr tidy --from-dir out/batch-2026-07/ -o analysis/
    """
    from frontstr.exports.tidy import (
        build_tidy_rows,
        load_payloads,
        parquet_available,
        write_tidy_csv,
        write_tidy_parquet,
    )

    paths = list(json_paths or [])
    if from_dir is not None:
        # Exclude the compact variant so a sample is not counted twice.
        paths += sorted(p for p in from_dir.rglob("*.json") if not p.name.endswith(".min.json"))
    if not paths:
        console.print("[red]error:[/red] no run JSONs given. Use --json or --from-dir.")
        raise typer.Exit(code=2)

    try:
        rows = build_tidy_rows(load_payloads(paths))
    except FrontstrError as exc:
        console.print(f"[red]tidy error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    samples = {r["sample"] for r in rows}
    out_dir.mkdir(parents=True, exist_ok=True)
    written = [write_tidy_csv(rows, out_dir / f"{stem}.csv")]
    if parquet_available():
        written.append(write_tidy_parquet(rows, out_dir / f"{stem}.parquet"))
    else:
        console.print(
            "[yellow]note:[/yellow] pyarrow not installed — CSV only. "
            "For the columnar dataset: pip install 'frontstr[parquet]'"
        )

    console.print(
        f"[green]{len(rows)} rows[/green] from {len(samples)} sample(s) across {len(paths)} run(s)"
    )
    for path in written:
        console.print(f"[green]wrote[/green] {path}  ({path.stat().st_size / 1024:.1f} KB)")


@app.command("calibrate-stutter")
def calibrate_stutter(
    bams: Annotated[
        list[Path],
        typer.Option("--bam", help="Indexed BAM/CRAM. Repeat for each sample."),
    ],
    panel_path: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    out: Annotated[
        Path | None, typer.Option("--out", "-o", help="Write the fitted model as JSON.")
    ] = None,
    protocol: Annotated[
        str,
        typer.Option(
            "--protocol",
            help="Library protocol the samples came from: wgs_pcr_free | amplicon. "
            "Recorded in the model — a PCR-free fit has no PCR slippage in it.",
        ),
    ] = "unknown",
    min_mapq: Annotated[int, typer.Option("--min-mapq")] = 20,
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", help="Reference FASTA (required for CRAM input)."),
    ] = None,
) -> None:
    """Measure stutter rates from real data and fit a :class:`StutterModel`.

    This is a calibration pass, not a per-case step: the model is a property of
    the chemistry and the library protocol, not of the sample. Run it once per
    platform + protocol, commit the JSON, and pass it to the pipeline.

    Only loci where a stutter position cannot be confused with a real allele
    are used (homozygotes, or heterozygotes ≥3 repeat units apart), so expect
    roughly half the loci to be discarded. Prints the breakdown by step, by LUS
    and by marker so a thin fit is visible rather than silently shipped.
    """
    from frontstr.panel.loader import load_panel
    from frontstr.panel.stutter_calib import (
        collect_observations,
        dump_stutter_model,
        fit_stutter_model,
        summarise,
    )

    try:
        panel = load_panel(panel_path)
    except FrontstrError as exc:
        console.print(f"[red]panel error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    obs = collect_observations(list(bams), panel, min_mapq=min_mapq, reference_fasta=reference)
    if not obs:
        console.print("[red]no usable observations[/red] — every locus was ambiguous.")
        raise typer.Exit(code=1)

    stats = summarise(obs)
    by_step: dict[str, dict[str, float]] = stats["by_step"]  # type: ignore[assignment]
    t = Table(title=f"Observed stutter — {len(bams)} sample(s), {panel.name}")
    t.add_column("Step", style="bold")
    t.add_column("n", justify="right")
    t.add_column("Pooled ratio", justify="right")
    t.add_column("Zero-stutter", justify="right")
    for step, row in by_step.items():
        t.add_row(step, str(row["n"]), f"{row['pooled_ratio']:.4f}", str(row["n_zero"]))
    console.print(t)

    lus_rows: dict[int, dict[str, float]] = stats["minus1_by_lus"]  # type: ignore[assignment]
    t2 = Table(title="-1 step by parent LUS (the covariate being fitted)")
    t2.add_column("LUS", justify="right")
    t2.add_column("n", justify="right")
    t2.add_column("Pooled ratio", justify="right")
    for lus, row in lus_rows.items():
        style = "" if row["n"] >= 7 else "dim"
        t2.add_row(str(lus), str(row["n"]), f"{row['pooled_ratio']:.4f}", style=style)
    console.print(t2)
    console.print("[dim]dim rows have thin support (n < 7) — widen the sample set.[/dim]")

    try:
        model = fit_stutter_model(
            obs,
            source=f"{len(bams)} sample(s): {', '.join(b.name for b in bams)}",
            protocol=protocol,
        )
    except ValueError as exc:
        console.print(f"[red]fit failed:[/red] {exc}")
        raise typer.Exit(code=1) from exc

    console.print(
        f"\n[green]rate(-1) = exp({model.log_intercept:+.4f} "
        f"+ {model.log_slope:.4f} * LUS)[/green]   "
        f"R²={model.r_squared}  LUS calibrado {model.lus_min}-{model.lus_max}"
    )
    console.print(f"step factors: {model.step_factors}   n={model.n_observations}")
    if protocol == "unknown":
        console.print(
            "[yellow]warning:[/yellow] --protocol not set. A model fitted on "
            "PCR-free data has no PCR slippage component and will under-predict "
            "stutter on amplicon casework."
        )

    if out is not None:
        dump_stutter_model(model, out)
        console.print(f"[green]✓ Model written:[/green] {out}")


@app.command("interpret")
def interpret(
    bam: Annotated[Path, typer.Option("--bam", help="Indexed sample BAM or CRAM.")],
    panel_path: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    longtr_vcf: Annotated[
        Path | None, typer.Option("--longtr-vcf", help="Optional LongTR VCF for concordance.")
    ] = None,
    min_mapq: Annotated[int, typer.Option("--min-mapq")] = 20,
    identity: Annotated[float, typer.Option("--identity")] = 0.97,
    len_tolerance: Annotated[int, typer.Option("--len-tolerance")] = 0,
    analytical_thresh: Annotated[
        float, typer.Option("--analytical-thresh", help="Fraction below = noise.")
    ] = 0.02,
    calling_thresh: Annotated[
        float, typer.Option("--calling-thresh", help="Fraction below = artefact.")
    ] = 0.10,
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", help="Reference FASTA (required for CRAM input)."),
    ] = None,
    catalog_path: Annotated[
        Path | None,
        typer.Option(
            "--catalog", help="Optional allele catalog JSON for ISFG/iso-allele annotation."
        ),
    ] = None,
) -> None:
    """End-to-end forensic call: pileup → cluster → ISFG → classify → call.

    This is the canonical FRONTStr command. It runs the evidence layer for
    every marker in ``--panel`` and prints one line per called allele plus
    LongTR concordance flags when ``--longtr-vcf`` is supplied.
    """
    from frontstr.caller import parse_longtr_vcf
    from frontstr.interp import index_longtr_results, interpret_run
    from frontstr.panel.catalog import load_catalog
    from frontstr.panel.loader import load_panel

    try:
        panel = load_panel(panel_path)
        catalog = load_catalog(catalog_path) if catalog_path else None
        longtr_map = index_longtr_results(parse_longtr_vcf(longtr_vcf)) if longtr_vcf else None
        results = interpret_run(
            bam=bam,
            panel=panel,
            longtr_results=longtr_map,
            min_mapq=min_mapq,
            identity_threshold=identity,
            len_tolerance_bp=len_tolerance,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
            reference_fasta=reference,
            catalog=catalog,
        )
    except FrontstrError as exc:
        console.print(f"[red]interpret error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    t = Table(title=f"FRONTStr — {len(results)} markers", show_lines=False)
    t.add_column("Marker", style="bold")
    t.add_column("Call")
    t.add_column("Tri", style="yellow")
    t.add_column("Alleles called")
    t.add_column("Cov", justify="right")
    t.add_column("LongTR?")
    for r in results:
        called = ", ".join(_interpret_allele_cell(a) for a in r.alleles_called)
        longtr_flag = (
            "[red]discordant[/red]" if r.discordant else ("ok" if r.longtr_result else "-")
        )
        t.add_row(
            r.marker_name,
            r.call_rule.value,
            r.tri_type.value or "-",
            called or "-",
            str(r.total_reads),
            longtr_flag,
        )
    console.print(t)


@app.command("evidence")
def evidence(
    bam: Annotated[Path, typer.Option("--bam", help="Indexed BAM or CRAM file.")],
    chrom: Annotated[str, typer.Option("--chrom", help="Chromosome (matches BAM @SQ).")],
    start: Annotated[int, typer.Option("--start", help="1-based inclusive TR start.")],
    end: Annotated[int, typer.Option("--end", help="1-based inclusive TR end.")],
    min_mapq: Annotated[int, typer.Option("--min-mapq", help="Drop reads below MAPQ.")] = 20,
    identity: Annotated[
        float, typer.Option("--identity", help="Cluster identity threshold (0-1).")
    ] = 0.97,
    len_tolerance: Annotated[
        int, typer.Option("--len-tolerance", help="bp tolerance for length binning.")
    ] = 0,
    motif: Annotated[
        str,
        typer.Option(
            "--motif",
            help="Comma-separated marker motifs (e.g. TCTA,TCTG). Enables "
            "repeat-core binning — pass it to reproduce what `interpret` does.",
        ),
    ] = "",
    strand: Annotated[
        str, typer.Option("--strand", help="Motif strand relative to the reference: + or -.")
    ] = "+",
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", help="Reference FASTA (required for CRAM input)."),
    ] = None,
) -> None:
    """Per-locus sequence pileup → clusters. Debug helper for the evidence layer.

    Prints one row per cluster (allele candidate) with integer read count and
    haplotype partition. This is the **forensic source of truth** for coverage
    in FRONTStr.
    """
    from frontstr.evidence.cluster import cluster_observations
    from frontstr.evidence.pileup import pileup_locus

    try:
        obs = pileup_locus(bam, chrom, start - 1, end, min_mapq=min_mapq, reference_fasta=reference)
    except FrontstrError as exc:
        console.print(f"[red]pileup error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    try:
        clusters = cluster_observations(
            obs,
            identity_threshold=identity,
            len_tolerance_bp=len_tolerance,
            motifs=[m for m in motif.split(",") if m],
            strand=strand,
        )
    except FrontstrError as exc:
        console.print(f"[red]cluster error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    title = f"Evidence: {chrom}:{start}-{end} ({len(obs)} obs, {len(clusters)} clusters)"
    t = Table(title=title)
    t.add_column("#", style="bold")
    t.add_column("Consensus")
    t.add_column("Len", justify="right")
    t.add_column("Reads", justify="right")
    t.add_column("HP1/HP2/none", justify="right")
    t.add_column("Strand +/-", justify="right")
    t.add_column("Mean Q", justify="right")
    for i, c in enumerate(clusters):
        cons = c.consensus
        cons_display = cons if len(cons) <= 60 else f"{cons[:28]}…{cons[-28:]}"
        t.add_row(
            str(i),
            cons_display,
            str(len(cons)),
            str(c.n_reads),
            f"{c.n_hp1}/{c.n_hp2}/{c.n_hp_none}",
            f"{c.n_forward}/{c.n_reverse}",
            f"{c.mean_qual:.1f}",
        )
    console.print(t)


if __name__ == "__main__":  # pragma: no cover
    app()
