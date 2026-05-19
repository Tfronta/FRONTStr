"""FRONTStr command-line interface.

Stages of the CLI are organised as Typer subcommands so each can be wired up
independently. Phase 1 ships the ``inspect`` subcommand (for debugging input
files) and a stub ``run`` that orchestrates the full pipeline.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from frontstr.errors import FrontstrError
from frontstr.ingest import detect_input, validate_bam
from frontstr.version import __version__

app = typer.Typer(
    name="frontstr",
    help="Forensic Ranked Output for Nanopore Tandem STR profiling.",
    no_args_is_help=True,
    rich_markup_mode="rich",
)
console = Console()


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"[bold]FRONTStr[/bold] {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    _version: Annotated[
        bool | None,
        typer.Option(
            "--version", "-V",
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
    panel: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    reference: Annotated[Path, typer.Option("--reference", "-r", help="Reference FASTA.")],
    out: Annotated[Path, typer.Option("--out", "-o", help="Output directory.")],
) -> None:
    """Run FRONTStr on a multi-sample batch (Phase 4)."""
    _ = manifest, panel, reference, out
    console.print(
        "[yellow]not implemented yet[/yellow] — Phase 4 of ROADMAP.md."
    )
    raise typer.Exit(code=64)


@app.command("call")
def call(
    bam: Annotated[Path, typer.Option("--bam", help="Indexed sample BAM.")],
    panel_path: Annotated[Path, typer.Option("--panel", "-p", help="Panel YAML.")],
    reference: Annotated[Path, typer.Option("--reference", "-r", help="Reference FASTA (indexed).")],
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
                panel=panel, reference=reference, platform=platform,
                phased=phased, binary=binary,
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
        gt = (
            "/".join(str(i) for i in sample.gt_indices)
            if sample and sample.gt_indices
            else "."
        )
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

    from frontstr.evidence.pileup import pileup_locus
    from frontstr.panel.loader import load_panel

    try:
        panel = load_panel(panel_path)
    except FrontstrError as exc:
        console.print(f"[red]panel error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    is_cram = bam.suffix.lower() == ".cram"
    if is_cram and reference is None:
        console.print("[red]error:[/red] CRAM input requires --reference <fasta>")
        raise typer.Exit(code=2)
    open_kwargs: dict = {"reference_filename": str(reference)} if is_cram else {}
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
            console.print(
                "[yellow]Hint:[/yellow] BAM has 'chr' prefix; panel does not."
            )

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
                    bam, s.chromosome, s.ref_start - 1, s.ref_end, min_mapq=0,
                    reference_fasta=reference,
                )
                obs_filt = pileup_locus(
                    bam, s.chromosome, s.ref_start - 1, s.ref_end, min_mapq=min_mapq,
                    reference_fasta=reference,
                )
                n_raw, n_filtered = len(obs_raw), len(obs_filt)
            except FrontstrError as exc:
                t.add_row(s.name, s.chromosome, f"{s.ref_start}-{s.ref_end}",
                          "yes", "?", "?", f"[red]{exc}[/red]")
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
            s.name, s.chromosome, f"{s.ref_start}-{s.ref_end}",
            "yes" if chrom_ok else "no",
            str(n_raw), str(n_filtered), status,
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
                "profile,evidence,seqs,json,json-compact,html."
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
) -> None:
    """Run the full pipeline and write one or more export files.

    Example::

        frontstr export --bam s.bam --panel p.yaml -o exports/ \
            --formats profile,evidence,seqs,json,html
    """
    from frontstr.caller import parse_longtr_vcf
    from frontstr.exports import (
        write_evidence_csv,
        write_profile_csv,
        write_run_json,
        write_seqs_csv,
    )
    from frontstr.interp import index_longtr_results, interpret_run
    from frontstr.panel.loader import load_panel
    from frontstr.report import RunContext, build_report, serialize_run

    wanted = {f.strip().lower() for f in formats.split(",") if f.strip()}
    unknown = wanted - {"profile", "evidence", "seqs", "json", "json-compact", "html"}
    if unknown:
        console.print(f"[red]Unknown formats:[/red] {sorted(unknown)}")
        raise typer.Exit(code=2)

    try:
        panel = load_panel(panel_path)
        longtr_map = (
            index_longtr_results(parse_longtr_vcf(longtr_vcf)) if longtr_vcf else None
        )
        results = interpret_run(
            bam=bam, panel=panel, longtr_results=longtr_map,
            min_mapq=min_mapq, identity_threshold=identity,
            analytical_thresh=analytical_thresh, calling_thresh=calling_thresh,
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
        payload = serialize_run(results, context)
    except FrontstrError as exc:
        console.print(f"[red]export error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    out_dir.mkdir(parents=True, exist_ok=True)
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
        written.append(
            write_run_json(payload, out_dir / f"{stem}.min.json", mode="compact")
        )
    if "html" in wanted:
        written.append(build_report(results, context, out_dir / f"{stem}.html"))

    for path in written:
        size_kb = path.stat().st_size / 1024
        console.print(f"[green]wrote[/green] {path}  ({size_kb:.1f} KB)")


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
        longtr_map = (
            index_longtr_results(parse_longtr_vcf(longtr_vcf)) if longtr_vcf else None
        )
        results = interpret_run(
            bam=bam, panel=panel, longtr_results=longtr_map,
            min_mapq=min_mapq, identity_threshold=identity,
            analytical_thresh=analytical_thresh, calling_thresh=calling_thresh,
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
) -> None:
    """End-to-end forensic call: pileup → cluster → ISFG → classify → call.

    This is the canonical FRONTStr command. It runs the evidence layer for
    every marker in ``--panel`` and prints one line per called allele plus
    LongTR concordance flags when ``--longtr-vcf`` is supplied.
    """
    from frontstr.caller import parse_longtr_vcf
    from frontstr.interp import index_longtr_results, interpret_run
    from frontstr.panel.loader import load_panel

    try:
        panel = load_panel(panel_path)
        longtr_map = (
            index_longtr_results(parse_longtr_vcf(longtr_vcf)) if longtr_vcf else None
        )
        results = interpret_run(
            bam=bam, panel=panel, longtr_results=longtr_map,
            min_mapq=min_mapq, identity_threshold=identity,
            len_tolerance_bp=len_tolerance,
            analytical_thresh=analytical_thresh, calling_thresh=calling_thresh,
            reference_fasta=reference,
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
        called = ", ".join(
            f"CE{a.ce}({a.n_reads_total})" if a.ce is not None else f"len{a.length_bp}({a.n_reads_total})"
            for a in r.alleles_called
        )
        longtr_flag = (
            "[red]discordant[/red]" if r.discordant
            else ("ok" if r.longtr_result else "-")
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
        obs = pileup_locus(bam, chrom, start - 1, end, min_mapq=min_mapq,
                           reference_fasta=reference)
    except FrontstrError as exc:
        console.print(f"[red]pileup error:[/red] {exc}")
        raise typer.Exit(code=2) from exc

    try:
        clusters = cluster_observations(
            obs, identity_threshold=identity, len_tolerance_bp=len_tolerance
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
