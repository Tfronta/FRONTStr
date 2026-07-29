"""Fetch panel-sized slices for a whole cohort out of the 1000G ONT bucket.

The published BAMs are whole-genome and tens of gigabytes each, so nothing here
downloads one. ``samtools`` is pointed at the remote URL together with its
remote index and asked for the panel regions only, which pulls the byte ranges
those regions occupy and nothing else — a few tens of megabytes per sample.

Two steps, kept separate on purpose so the expensive one runs once::

    python -m validation.cohort fetch --out cohort/slices
    frontstr batch --manifest cohort/slices/manifest.tsv \\
        -p examples/panels/codis_20_grch38.yaml -o cohort/run
    python -m validation.compare --summary cohort/run/batch_summary.csv \\
        --workbook ~/Desktop/1000GEN-ONT-Merged-Compar.xlsx --out cohort

``fetch`` is resumable: a sample whose slice and index already exist is skipped,
so an interrupted run continues where it stopped. Keeping the slices means a
caller change can be re-scored across the cohort in minutes instead of hours.

Only R10 + Dorado basecalls are eligible. R9 and guppy are rejected outright —
the stutter model and every ``corr_value`` in the panel were calibrated on R10 +
Dorado, so an R9 sample would be scored against a model that does not describe
it.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from frontstr.panel.loader import load_panel
from validation.truth import canonical_sample_id, load_truth

app = typer.Typer(add_completion=False, help=__doc__)
console = Console()

#: Public HTTPS front door for ``s3://1000g-ont``. Anonymous, no credentials.
BUCKET_URL = "https://s3.amazonaws.com/1000g-ont"

#: Only these basecalls may be used. See the ``ont-r10-dorado-only`` rule.
REQUIRED_IN_FILENAME = ("-R10", "dorado")
REJECTED_IN_FILENAME = ("-R9", "guppy")

#: Default bucket path listing: one remote BAM path per line.
DEFAULT_PATH_LIST = Path.home() / "hg38_minimap2_align_paths.txt"

#: Padding around each panel window, in bases. A read must span the repeat plus
#: its flanks to be usable, so the slice has to be wider than the window itself.
DEFAULT_PADDING = 10_000


@dataclass(frozen=True, slots=True)
class RemoteSample:
    """One cohort sample and the remote BAM it will be sliced from."""

    sample_id: str
    remote_path: str

    @property
    def url(self) -> str:
        return f"{BUCKET_URL}/{self.remote_path}"

    @property
    def index_url(self) -> str:
        return f"{self.url}.bai"


def eligible_basecall(remote_path: str) -> bool:
    """True when the filename says R10 + Dorado and neither R9 nor guppy."""
    name = remote_path.rsplit("/", 1)[-1]
    if any(bad.lower() in name.lower() for bad in REJECTED_IN_FILENAME):
        return False
    return all(good.lower() in name.lower() for good in REQUIRED_IN_FILENAME)


def parse_path_list(path: Path) -> dict[str, list[str]]:
    """Group the bucket path listing by sample ID, keeping only R10 + Dorado.

    The sample ID is the filename up to the first ``-ONT``; ``GM`` IDs are
    rewritten to the ``NA`` spelling the workbook uses.
    """
    if not path.exists():
        raise FileNotFoundError(f"Path listing not found: {path}")

    by_sample: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        remote = line.strip()
        if not remote or remote.startswith("#") or not eligible_basecall(remote):
            continue
        name = remote.rsplit("/", 1)[-1]
        sample_id = canonical_sample_id(name.split("-ONT", 1)[0])
        by_sample.setdefault(sample_id, []).append(remote)
    return by_sample


def select_cohort(
    path_list: Path,
    workbook: Path | None,
) -> tuple[list[RemoteSample], list[str]]:
    """Pick the samples to run and report the ones with no eligible BAM.

    When ``workbook`` is given the cohort is exactly the samples that have truth
    genotypes — running the rest produces calls nobody can score.
    """
    available = parse_path_list(path_list)

    if workbook is not None:
        wanted = sorted({call.sample_id for call in load_truth(workbook)})
    else:
        wanted = sorted(available)

    selected: list[RemoteSample] = []
    missing: list[str] = []
    for sample_id in wanted:
        paths = sorted(available.get(sample_id, []))
        if not paths:
            missing.append(sample_id)
            continue
        # More than one eligible basecall is possible (different Dorado
        # versions); take the last sorted so the newest wins, deterministically.
        selected.append(RemoteSample(sample_id=sample_id, remote_path=paths[-1]))
    return selected, missing


def write_slice_bed(panel_path: Path, out_path: Path, padding: int) -> Path:
    """Write the padded slicing BED, derived from the panel itself.

    Derived rather than committed so the slice regions cannot drift away from
    the windows the caller will actually read.
    """
    panel = load_panel(panel_path)
    lines: list[str] = []
    for system in panel.systems:
        start = max(0, system.ref_start - 1 - padding)
        end = system.ref_end + padding
        lines.append(f"{system.chromosome}\t{start}\t{end}\t{system.name}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


def fetch_slice(sample: RemoteSample, bed: Path, out_dir: Path, *, timeout: int) -> Path:
    """Download just the panel regions of one remote BAM, and index them.

    ``-M`` makes ``-L`` a true multi-region filter rather than a union of
    overlapping fetches, and ``-X`` lets the remote index be named separately so
    htslib never scans the whole object.
    """
    out_bam = out_dir / f"{sample.sample_id}.codis.bam"
    if out_bam.exists() and (out_bam.with_suffix(".bam.bai")).exists():
        return out_bam

    partial = out_bam.with_suffix(".bam.partial")
    subprocess.run(
        [
            "samtools",
            "view",
            "-b",
            "-M",
            "-L",
            str(bed),
            "-X",
            sample.url,
            sample.index_url,
            "-o",
            str(partial),
        ],
        check=True,
        capture_output=True,
        timeout=timeout,
    )
    # Only claim the final name once the download completed, so an interrupted
    # run leaves no truncated BAM that the resume logic would happily skip.
    partial.replace(out_bam)
    subprocess.run(
        ["samtools", "index", str(out_bam)], check=True, capture_output=True, timeout=timeout
    )
    return out_bam


def write_manifest(samples: list[tuple[str, Path]], out_path: Path) -> Path:
    """Write the batch manifest ``frontstr batch --manifest`` consumes."""
    lines = ["# FRONTStr cohort validation manifest", "sample_id\tbam\trole"]
    lines += [f"{sample_id}\t{bam}\tsample" for sample_id, bam in samples]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


@app.command()
def fetch(
    out: Annotated[Path, typer.Option("--out", "-o", help="Directory for the slices")],
    path_list: Annotated[
        Path, typer.Option("--paths", help="Bucket path listing, one BAM per line")
    ] = DEFAULT_PATH_LIST,
    workbook: Annotated[
        Path | None, typer.Option("--workbook", "-w", help="Restrict to samples with truth")
    ] = None,
    panel: Annotated[
        Path, typer.Option("--panel", "-p", help="Panel YAML the BED is derived from")
    ] = Path("examples/panels/codis_20_grch38.yaml"),
    padding: Annotated[int, typer.Option("--padding", help="Bases around each window")] = (
        DEFAULT_PADDING
    ),
    limit: Annotated[int, typer.Option("--limit", help="Stop after N samples (0 = all)")] = 0,
    timeout: Annotated[int, typer.Option("--timeout", help="Per-sample seconds")] = 1800,
) -> None:
    """Download panel slices for the cohort and write the batch manifest."""
    out.mkdir(parents=True, exist_ok=True)
    bed = write_slice_bed(panel, out / "slice_regions.bed", padding)

    selected, missing = select_cohort(path_list, workbook)
    if limit:
        selected = selected[:limit]

    console.print(f"Cohort: [bold]{len(selected)}[/bold] sample(s), slicing with {bed}")
    if missing:
        console.print(
            f"[yellow]No R10+Dorado BAM for {len(missing)} sample(s):[/yellow] "
            f"{', '.join(missing[:10])}{' …' if len(missing) > 10 else ''}"
        )

    done: list[tuple[str, Path]] = []
    failed: list[tuple[str, str]] = []
    for n, sample in enumerate(selected, start=1):
        cached = (out / f"{sample.sample_id}.codis.bam").exists()
        console.print(
            f"  [{n}/{len(selected)}] {sample.sample_id}" + ("  (cached)" if cached else "")
        )
        try:
            bam = fetch_slice(sample, bed, out, timeout=timeout)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            detail = getattr(exc, "stderr", b"") or b""
            failed.append((sample.sample_id, detail.decode("utf-8", "replace").strip()[:200]))
            console.print(f"    [red]failed:[/red] {failed[-1][1] or type(exc).__name__}")
            continue
        done.append((sample.sample_id, bam.resolve()))

    manifest = write_manifest(done, out / "manifest.tsv")
    total_mb = sum(bam.stat().st_size for _, bam in done) / 1e6
    console.print(
        f"\n[bold]{len(done)}/{len(selected)}[/bold] slices ready ({total_mb:.0f} MB)"
        f"{f', {len(failed)} failed' if failed else ''}"
    )
    console.print(f"Manifest: [bold]{manifest}[/bold]")


@app.command()
def plan(
    path_list: Annotated[
        Path, typer.Option("--paths", help="Bucket path listing, one BAM per line")
    ] = DEFAULT_PATH_LIST,
    workbook: Annotated[
        Path | None, typer.Option("--workbook", "-w", help="Restrict to samples with truth")
    ] = None,
) -> None:
    """Report what ``fetch`` would download, without downloading anything."""
    selected, missing = select_cohort(path_list, workbook)
    console.print(f"Eligible (R10 + Dorado): [bold]{len(selected)}[/bold] sample(s)")
    if missing:
        console.print(f"[yellow]No eligible BAM:[/yellow] {len(missing)} — {', '.join(missing)}")
    for sample in selected[:5]:
        console.print(f"  {sample.sample_id}  {sample.remote_path.rsplit('/', 1)[-1]}")
    if len(selected) > 5:
        console.print(f"  … and {len(selected) - 5} more")


if __name__ == "__main__":
    app()
