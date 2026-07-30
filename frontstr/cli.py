"""FRONTStr command-line interface.

Stages of the CLI are organised as Typer subcommands so each can be wired up
independently. Phase 1 ships the ``inspect`` subcommand (for debugging input
files) and a stub ``run`` that orchestrates the full pipeline.
"""

from __future__ import annotations

import platform
import sys
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


def _load_regions(
    panel_path: Path | None, bed: Path | None, bed_coords: str
) -> tuple[Any, list[str]]:
    """Resolve ``--panel`` / ``--bed`` into a panel, plus warnings to print.

    A BED is the escape hatch HipSTR and LongTR have and a YAML-only caller does
    not: point the tool at your own regions without curating a panel first. What
    it cannot carry is calibration, so the warnings are not decoration — for a
    marker STRNaming has no range for, the number is an uncalibrated repeat
    count rather than a kit allele, and the run says so per locus as well.
    """
    from frontstr.panel.bed import load_panel_from_bed
    from frontstr.panel.loader import load_panel

    if bed is not None and panel_path is not None:
        raise FrontstrError("give --panel or --bed, not both")
    if bed is None and panel_path is None:
        raise FrontstrError("need --panel <yaml> or --bed <file>")
    if bed is None:
        assert panel_path is not None
        return load_panel(panel_path), []

    if bed_coords not in ("bed0", "panel1"):
        raise FrontstrError(f"--bed-coords must be bed0 or panel1, got {bed_coords!r}")
    panel, uncalibrated = load_panel_from_bed(bed, coords=bed_coords)  # type: ignore[arg-type]
    warnings = [
        f"[dim]Regions from {bed} read as "
        f"{'0-based half-open (standard BED)' if bed_coords == 'bed0' else '1-based inclusive'}"
        f" — {len(panel.systems)} marker(s).[/dim]"
    ]
    if uncalibrated:
        warnings.append(
            f"[yellow]warning:[/yellow] {len(uncalibrated)} marker(s) have no STRNaming "
            "reporting range and no calibration from the BED, so their allele number is "
            "an uncalibrated repeat count, not a kit CE allele: "
            f"{', '.join(uncalibrated[:8])}" + (" …" if len(uncalibrated) > 8 else "")
        )
    return panel, warnings


def _environment_report() -> list[str]:
    """Check the installation itself. Returns the problems found, if any.

    The failure this exists for is quiet: without a POA backend FRONTStr still
    emits a complete profile, built from unpolished single reads, and the damage
    surfaces as microvariants that are not in the sample. Measured on the
    reference slices, that fallback produced 4 false microvariants in 202 called
    alleles. Nobody should have to discover it from a flag after the fact.
    """
    from frontstr.evidence.consensus import poa_backend_name
    from frontstr.interp.naming import default_namer
    from frontstr.version import __version__

    problems: list[str] = []
    t = Table(title="FRONTStr doctor — environment", show_header=False)
    t.add_column("Check", style="bold")
    t.add_column("Result")

    t.add_row("FRONTStr", __version__)
    t.add_row("Python", sys.version.split()[0])
    t.add_row("Platform", f"{platform.system()} {platform.machine()}")

    backend = poa_backend_name()
    if backend:
        t.add_row("POA backend", f"[green]{backend}[/green]")
    else:
        t.add_row("POA backend", "[red]none — consensus will be a single unpolished read[/red]")
        problems.append("no POA backend: pip install 'frontstr[poa]'")

    namer = default_namer()
    if namer is None:
        t.add_row(
            "STRNaming", "[red]unavailable — allele numbers fall back to bracket counts[/red]"
        )
        problems.append("STRNaming unavailable: pip install 'strnaming>=1.2,<1.3'")
    else:
        n = sum(1 for _ in _cached_markers())
        t.add_row("STRNaming", f"[green]ready[/green], {n} markers in the bundled slice cache")

    for mod in ("pysam", "edlib", "cyvcf2"):
        try:
            m = __import__(mod)
            t.add_row(mod, getattr(m, "__version__", "installed"))
        except ImportError:
            t.add_row(mod, "[red]missing[/red]")
            problems.append(f"{mod} is not importable")

    console.print(t)
    for p in problems:
        console.print(f"  [red]✗[/red] {p}")
    if not problems:
        console.print("  [green]✓[/green] installation looks complete")
    return problems


def _cached_markers() -> list[str]:
    from frontstr.interp.naming import CACHE_PATH

    try:
        lines = CACHE_PATH.read_text(encoding="utf-8").splitlines()[1:]
    except OSError:
        return []
    return [x.split("\t")[0] for x in lines if x.strip()]


def _phasing_report(af: Any, panel: Any, reference: Path | None) -> None:
    """Whether the BAM is phased, and whether its blocks survive a locus.

    Haplotype evidence drives two calling rules, and both are silently disabled
    on an unphased BAM. It is better to know that before reading a profile than
    to wonder why a rescue never fired.
    """
    hp = ps = total = 0
    blocks_per_locus: list[int] = []
    for s in panel.systems[:8]:
        seen: set[int] = set()
        try:
            reads = af.fetch(s.chromosome, max(0, s.ref_start - 1), s.ref_end)
        except (ValueError, OSError):
            continue
        for r in reads:
            total += 1
            if r.has_tag("HP"):
                hp += 1
                if r.has_tag("PS"):
                    ps += 1
                    seen.add(int(r.get_tag("PS")))
        if seen:
            blocks_per_locus.append(len(seen))

    if total == 0:
        return
    if hp == 0:
        console.print(
            "[yellow]Phasing:[/yellow] no HP tags — haplotype rules are disabled "
            "(phantom suppression and the peak-ratio rescue both no-op)."
        )
        return
    pct = 100 * hp / total
    line = f"[green]Phasing:[/green] {pct:.0f}% of sampled reads carry HP"
    if ps == 0:
        line += "; [yellow]no PS tags[/yellow] — blocks cannot be verified"
    else:
        split = sum(1 for n in blocks_per_locus if n > 1)
        line += f", {ps} with PS"
        if split:
            line += f"; [yellow]{split} sampled locus/loci span >1 phase block[/yellow]"
    console.print(line)


def _render_params(params: Any) -> str:
    """The parameter block, with overridden values highlighted for a terminal."""
    from rich.markup import escape

    from frontstr.params import render_echo

    out = []
    for line in render_echo(params).splitlines():
        if "CHANGED" in line:
            out.append(f"[yellow]{escape(line)}[/yellow]")
        elif line.startswith("Parameters in force"):
            out.append(f"[bold]{escape(line)}[/bold]")
        else:
            out.append(f"[dim]{escape(line)}[/dim]")
    return "\n".join(out)


def _coverage_cell(result: Any) -> str:
    """Depth behind the call: the per-allele read counts, added up. Nothing else.

    This used to trail ``+n``, the spanning reads that supported no called
    allele. They were dropped: the caller had already discarded them, so
    carrying them in the depth column reports discarded evidence as depth. The
    same reasoning removed them from the HTML report, and the two views have to
    agree about what a number means.

    Where those reads went is not lost, it moved to where it can be acted on.
    ``--trace`` names every one of them and the rule that discarded it, which
    is the view for a locus where most reads support nothing.
    """
    return str(result.called_reads) if result.alleles_called else "0"


def _balance_cell(result: Any) -> str:
    """Allele balance, dimmed while inside the balanced band."""
    ab = result.allele_balance
    if ab is None:
        return "[dim]-[/dim]"
    thr = QcThresholds().balanced_ab_max
    return f"{ab:.2f}" if ab > thr else f"[dim]{ab:.2f}[/dim]"


def _qc_cell(result: Any) -> str:
    """The flags that actually fired, worst-coloured. Never an aggregated PASS.

    A green PASS standing for several checks teaches a reviewer to stop reading
    the individual ones, and one that shows on almost every locus stops carrying
    information. A clean locus therefore says nothing at all.
    """
    if not result.flags:
        return "[dim]-[/dim]"
    colour = {"error": "red", "warn": "yellow", "info": "cyan"}
    return " ".join(
        f"[{colour.get(f.severity.value, 'white')}]{f.code.value}[/]" for f in result.flags
    )


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
    layers (ingest → align/passthrough → evidence → interp → report).
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
    log: Annotated[
        bool, typer.Option("--log", help="Watch the run: one line per marker, per sample.")
    ] = False,
    trace: Annotated[
        bool,
        typer.Option(
            "--trace/--no-trace",
            help="Per-locus narrative, one file per sample. On by default.",
        ),
    ] = True,
) -> None:
    """Run FRONTStr on a multi-sample batch from a manifest TSV.

    Pass ``--log`` to follow the run rather than only its progress counter:
    the per-marker process log on stderr, every line tagged with the sample it
    came from, which is what makes it readable when ``-j`` runs several at once.

    Every run writes the full account of every locus to
    ``<out>/<sample>/<sample>.trace.txt`` — the read funnel with each rejection
    reason, the bins, the clusters with their consensus, the aligned sequences,
    the HP1/HP2 haplotype counts, how each allele was named, and the call.

    This is on by default and should stay on. It is the only record that lets
    someone ask *where* a call went wrong instead of only whether it did, and
    it costs nothing measurable: 18.5 s against 19.0 s over five samples,
    ~160 kB per sample. ``--no-trace`` exists for the case where the output
    directory is genuinely constrained; a run without it cannot be questioned
    later without being repeated.

    **With ``-j 1`` the narrative also appears on screen as it happens.** With
    more than one worker it does not — the loci of different samples would land
    interleaved and become unreadable — so there ``--log`` is the live channel
    and the trace is read afterwards.

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
            log=log,
            trace=trace,
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


@app.command("doctor")
def doctor(
    bam: Annotated[Path | None, typer.Option("--bam", help="Indexed sample BAM or CRAM.")] = None,
    panel_path: Annotated[Path | None, typer.Option("--panel", "-p", help="Panel YAML.")] = None,
    min_mapq: Annotated[int, typer.Option("--min-mapq")] = 20,
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", help="Reference FASTA (required for CRAM input)."),
    ] = None,
) -> None:
    """Pre-flight check: the environment, then BAM ↔ panel compatibility.

    With no arguments it checks the installation alone — the POA backend, the
    STRNaming slice cache, the compiled dependencies. Worth running after any
    install, because the failure mode is quiet: without a POA backend FRONTStr
    still produces a full profile, from unpolished single reads, and the damage
    shows up as microvariants that are not there.

    Given ``--bam`` and ``--panel`` it also runs every marker through the
    pileup: whether each chromosome exists in the BAM headers, read counts with
    and without the MAPQ filter, and whether the reads are phased.

    Exits non-zero if anything is broken, so it can gate a batch.
    """
    import pysam

    from frontstr.evidence.pileup import pileup_locus
    from frontstr.panel.loader import load_panel

    problems = _environment_report()

    if bam is None or panel_path is None:
        if bam is not None or panel_path is not None:
            console.print(
                "[yellow]note:[/yellow] --bam and --panel are checked together; "
                "give both for the per-marker table."
            )
        if problems:
            raise typer.Exit(code=1)
        return

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
    _phasing_report(af, panel, reference)
    af.close()
    if problems:
        raise typer.Exit(code=1)


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
    import sys

    from frontstr.exports import (
        write_evidence_csv,
        write_profile_csv,
        write_run_json,
        write_run_vcf,
        write_run_xlsx,
        write_seqs_csv,
    )
    from frontstr.interp import interpret_run
    from frontstr.log import PROCESS_LOG_NAME, configure_logging
    from frontstr.panel.bed import panel_bed_lines
    from frontstr.panel.loader import load_panel
    from frontstr.params import RunParameters
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
        results = interpret_run(
            bam=bam,
            panel=panel,
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
            platform=platform,
            operator=operator,
            run_id=run_id,
            reference_build=panel.reference_build,
            qc_thresholds=qc_thresholds,
            pipeline_argv=list(sys.argv),
            effective_params=RunParameters.of(
                min_mapq=min_mapq,
                identity_threshold=identity,
                analytical_thresh=analytical_thresh,
                calling_thresh=calling_thresh,
                low_coverage_reads=qc_thresholds.low_coverage_reads,
                balanced_ab_max=qc_thresholds.balanced_ab_max,
            ).as_audit_rows(),
            panel_bed=panel_bed_lines(panel),
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

        frontstr report --bam s.bam --panel codis.yaml --out s.html
    """
    import sys

    from frontstr.interp import interpret_run
    from frontstr.panel.bed import panel_bed_lines
    from frontstr.panel.loader import load_panel
    from frontstr.params import RunParameters
    from frontstr.report import RunContext, build_report

    try:
        panel = load_panel(panel_path)
        results = interpret_run(
            bam=bam,
            panel=panel,
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
            platform=platform,
            operator=operator,
            run_id=run_id,
            reference_build=panel.reference_build,
            pipeline_argv=list(sys.argv),
            effective_params=RunParameters.of(
                min_mapq=min_mapq,
                identity_threshold=identity,
                analytical_thresh=analytical_thresh,
                calling_thresh=calling_thresh,
            ).as_audit_rows(),
            panel_bed=panel_bed_lines(panel),
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
    panel_path: Annotated[Path | None, typer.Option("--panel", "-p", help="Panel YAML.")] = None,
    bed: Annotated[
        Path | None,
        typer.Option("--bed", help="Regions as BED instead of a panel. See --bed-coords."),
    ] = None,
    bed_coords: Annotated[
        str,
        typer.Option(
            "--bed-coords",
            help="How to read --bed: bed0 (standard 0-based half-open) or panel1 (1-based inclusive).",
        ),
    ] = "bed0",
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
    flank_anchor: Annotated[
        int, typer.Option("--flank-anchor", help="bp of clean flank required each side.")
    ] = 20,
    min_phr_for_het: Annotated[
        float,
        typer.Option("--min-phr", help="Minor/major read ratio to call a heterozygote."),
    ] = 0.4,
    min_reads_third: Annotated[
        int | None,
        typer.Option(
            "--min-reads-third", help="Read floor for a 3rd allele. DERIVED — see --help."
        ),
    ] = None,
    low_coverage_reads: Annotated[
        int,
        typer.Option("--low-coverage-reads", help="Flag a call below this coverage. DERIVED."),
    ] = 20,
    balanced_ab_max: Annotated[
        float, typer.Option("--balanced-ab-max", help="Largest balanced allele-balance value.")
    ] = 0.65,
    log: Annotated[
        bool,
        typer.Option("--log", "-l", help="Print the per-marker process log to stderr."),
    ] = False,
    show_params: Annotated[
        bool,
        typer.Option(
            "--show-params/--no-show-params",
            help="Print every parameter in force before the run.",
        ),
    ] = True,
    trace: Annotated[
        bool,
        typer.Option("--trace", help="Narrate every pipeline step, per locus, to stderr."),
    ] = False,
    trace_out: Annotated[
        Path | None,
        typer.Option("--trace-out", help="Write the narrative trace to a file instead."),
    ] = None,
) -> None:
    """End-to-end forensic call: pileup → cluster → ISFG → classify → call.

    This is the canonical FRONTStr command. It runs the evidence layer for
    every marker in ``--panel`` and prints one line per called allele.

    Pass ``--log`` to watch what it does rather than only what it concluded:
    the run configuration, then one line per marker with the call rule, the
    coverage, the cluster count and how each allele number was derived.

    Pass ``--trace`` for the full account of every locus: the read funnel with
    each rejection reason, the repeat-core bins, the clusters and their
    consensus, how each candidate was named and classified, and the call. This
    is the transparent-validation view — a genotype can be followed all the way
    back to the reads.

    Both go to stderr, so ``frontstr interpret ... --trace > table.txt`` still
    separates cleanly. ``--trace-out FILE`` writes the narrative to a file.
    """
    import logging
    import sys

    from frontstr.interp import interpret_run
    from frontstr.log import configure_logging
    from frontstr.panel.catalog import load_catalog
    from frontstr.params import RunParameters
    from frontstr.trace import (
        LocusTrace,
        RunHeader,
        render_header,
        render_locus,
        render_run_summary,
    )

    if log:
        configure_logging(level=logging.DEBUG, console=True)

    params = RunParameters.of(
        min_mapq=min_mapq,
        flank_anchor=flank_anchor,
        identity_threshold=identity,
        len_tolerance_bp=len_tolerance,
        analytical_thresh=analytical_thresh,
        calling_thresh=calling_thresh,
        min_phr_for_het=min_phr_for_het,
        min_reads_third=min_reads_third,
        low_coverage_reads=low_coverage_reads,
        balanced_ab_max=balanced_ab_max,
    )
    if show_params:
        console.print(_render_params(params))

    traces: list[LocusTrace] = []
    trace_fh = trace_out.open("w", encoding="utf-8") if trace_out else None
    want_trace = trace or trace_out is not None

    def emit(locus: LocusTrace) -> None:
        traces.append(locus)
        text = render_locus(locus)
        if trace_fh is not None:
            trace_fh.write(text + "\n\n")
        else:
            print(text, file=sys.stderr)
            print(file=sys.stderr)

    try:
        panel, bed_warnings = _load_regions(panel_path, bed, bed_coords)
        for warning in bed_warnings:
            console.print(warning)
        catalog = load_catalog(catalog_path) if catalog_path else None
        if want_trace:
            from frontstr.evidence.consensus import poa_backend_name
            from frontstr.interp.naming import default_namer
            from frontstr.panel.stutter_calib import DEFAULT_STUTTER_MODEL
            from frontstr.version import __version__ as _v

            namer = default_namer()
            head = render_header(
                RunHeader(
                    inputs=[str(bam)],
                    panel_name=panel.name,
                    panel_version=panel.version or "",
                    n_markers=len(panel.systems),
                    min_mapq=min_mapq,
                    identity_threshold=identity,
                    analytical_thresh=analytical_thresh,
                    calling_thresh=calling_thresh,
                    consensus_backend=poa_backend_name(),
                    naming_markers=sum(
                        1 for s in panel.systems if namer and namer.has_range(s.name)
                    ),
                    stutter_model=DEFAULT_STUTTER_MODEL.describe(),
                    tool_version=_v,
                    overrides=[
                        (spec.name, params[spec.name], spec.default, spec.provenance)
                        for spec in params.overrides()
                    ],
                )
            )
            if trace_fh is not None:
                trace_fh.write(head + "\n")
            else:
                print(head, file=sys.stderr)
        results = interpret_run(
            bam=bam,
            panel=panel,
            reference_fasta=reference,
            catalog=catalog,
            on_trace=emit if want_trace else None,
            params=params,
            qc_thresholds=QcThresholds(
                low_coverage_reads=low_coverage_reads, balanced_ab_max=balanced_ab_max
            ),
        )
    except FrontstrError as exc:
        console.print(f"[red]interpret error:[/red] {exc}")
        raise typer.Exit(code=2) from exc
    finally:
        if want_trace and traces:
            summary = render_run_summary(traces)
            if trace_fh is not None:
                trace_fh.write(summary + "\n")
            else:
                print(summary, file=sys.stderr)
        if trace_fh is not None:
            trace_fh.close()
            console.print(f"[dim]trace written to {trace_out}[/dim]")

    t = Table(title=f"FRONTStr — {len(results)} markers", show_lines=False)
    t.add_column("Marker", style="bold")
    t.add_column("Call")
    t.add_column("Tri", style="yellow")
    # no_wrap: an allele list broken across two lines costs more room than the
    # column saves, and makes the per-allele reads hard to pair with a number.
    t.add_column("Alleles called", no_wrap=True)
    t.add_column("Cov", justify="right")
    t.add_column("AB", justify="right")
    # fold, not truncate: a QC column that shows "allele_imbala…" is worse
    # than one that wraps, because the reader cannot tell which flag fired.
    t.add_column("QC", overflow="fold")
    for r in results:
        called = ", ".join(_interpret_allele_cell(a) for a in r.alleles_called)
        t.add_row(
            r.marker_name,
            r.call_rule.value,
            r.tri_type.value or "-",
            called or "-",
            _coverage_cell(r),
            _balance_cell(r),
            _qc_cell(r),
        )
    console.print(t)
    console.print(
        "[dim]Alleles called: number(reads on that allele). Cov: the depth behind the "
        "call, those per-allele counts added up. Reads that supported no called allele "
        "were discarded by the caller and are not counted here; run --trace to see what "
        "they were and which rule discarded each.[/dim]"
    )


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
