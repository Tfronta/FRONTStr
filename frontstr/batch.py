"""Batch multi-sample orchestrator.

Usage
-----
High-level entry point for running FRONTStr on a cohort of samples described
by a tab-separated manifest file::

    sample_id   bam                      role
    HG00113     /data/HG00113.bam        sample
    CTRL_POS    /data/positive.bam       positive_ctrl
    CTRL_NEG    /data/negative.bam       negative_ctrl
    BLANK_1     /data/blank.bam          reagent_blank

Manifest rules
--------------
- First non-comment line must be the header: ``sample_id``, ``bam``, and
  optionally ``role``.
- Lines starting with ``#`` are ignored.
- ``role`` defaults to ``sample`` when the column is absent or empty.
- Valid roles: ``sample``, ``positive_ctrl``, ``negative_ctrl``,
  ``reagent_blank``.

Output layout
-------------
::

    out_dir/
      HG00113/
        HG00113.profile.csv
        HG00113.evidence.csv
        HG00113.seqs.csv
        HG00113.json
        HG00113.html
        HG00113.trace.txt      <- per-locus narrative; see frontstr.audit
      CTRL_POS/
        ...
      batch_summary.csv
"""

from __future__ import annotations

import csv
import logging
import sys
import traceback
from collections.abc import Callable, Iterator
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from frontstr.errors import FrontstrError
from frontstr.interp.models import MarkerResult
from frontstr.panel.models import Panel

VALID_ROLES: frozenset[str] = frozenset(
    {"sample", "positive_ctrl", "negative_ctrl", "reagent_blank"}
)

DEFAULT_FORMATS: frozenset[str] = frozenset({"profile", "evidence", "seqs", "json", "html"})

SUMMARY_BASE_HEADERS = ("sample_id", "role", "status", "error")


@dataclass(slots=True)
class ManifestEntry:
    """One row from the batch manifest."""

    sample_id: str
    bam: Path
    role: str = "sample"


@dataclass(slots=True)
class BatchResult:
    """Outcome for one sample in a batch run."""

    sample_id: str
    role: str
    status: str  # "ok" | "error"
    error: str = ""
    files: list[Path] = field(default_factory=list)
    marker_ces: dict[str, str] = field(default_factory=dict)


def parse_manifest(path: Path) -> list[ManifestEntry]:
    """Parse a tab-separated batch manifest and return :class:`ManifestEntry` list.

    Args:
        path: Path to the manifest TSV.

    Raises:
        FrontstrError: On missing file, bad header, duplicate sample IDs, or
            unknown role values.
    """
    if not path.exists():
        raise FrontstrError(f"Manifest not found: {path}")

    entries: list[ManifestEntry] = []
    seen_ids: set[str] = set()

    with path.open(encoding="utf-8", newline="") as fh:
        # Skip comment lines and find the header
        lines = [ln.rstrip("\n") for ln in fh]

    data_lines = [ln for ln in lines if ln and not ln.startswith("#")]
    if not data_lines:
        raise FrontstrError(f"Manifest {path} is empty (no non-comment lines)")

    reader = csv.DictReader(data_lines, delimiter="\t")
    if (
        reader.fieldnames is None
        or "sample_id" not in reader.fieldnames
        or "bam" not in reader.fieldnames
    ):
        raise FrontstrError(
            f"Manifest {path} must have tab-separated columns: sample_id, bam[, role]"
        )

    for lineno, row in enumerate(reader, start=2):
        sid = row.get("sample_id", "").strip()
        bam_str = row.get("bam", "").strip()
        role = (row.get("role") or "").strip() or "sample"

        if not sid:
            raise FrontstrError(f"Manifest line {lineno}: empty sample_id")
        if not bam_str:
            raise FrontstrError(f"Manifest line {lineno}: empty bam path")
        if role not in VALID_ROLES:
            raise FrontstrError(
                f"Manifest line {lineno}: unknown role {role!r}; "
                f"valid roles are {sorted(VALID_ROLES)}"
            )
        if sid in seen_ids:
            raise FrontstrError(f"Manifest line {lineno}: duplicate sample_id {sid!r}")

        seen_ids.add(sid)
        entries.append(ManifestEntry(sample_id=sid, bam=Path(bam_str), role=role))

    if not entries:
        raise FrontstrError(f"Manifest {path} has a header but no data rows")

    return entries


def run_batch(
    *,
    entries: list[ManifestEntry],
    panel: Panel,
    out_dir: Path,
    reference_fasta: Path | None = None,
    formats: frozenset[str] = DEFAULT_FORMATS,
    min_mapq: int = 20,
    identity: float = 0.97,
    analytical_thresh: float = 0.02,
    calling_thresh: float = 0.10,
    platform: str = "ont",
    operator: str | None = None,
    run_id: str | None = None,
    workers: int = 1,
    progress_callback: Any | None = None,
    log: bool = False,
    trace: bool = True,
) -> list[BatchResult]:
    """Run the full pipeline on each manifest entry, optionally in parallel.

    Args:
        entries: Parsed manifest entries.
        panel: Loaded panel definition.
        out_dir: Root output directory; per-sample subdirectories are created.
        reference_fasta: Reference FASTA path; required for CRAM entries.
        formats: Set of output formats to write per sample.
        workers: Number of parallel processes (1 = sequential, no subprocess).
        progress_callback: Optional callable invoked as ``callback(sample_id)``
            after each sample completes (used by CLI for live progress).
        log: Emit the per-marker process log on stderr as each sample runs,
            every line tagged with its sample.
        trace: Write the full per-locus narrative to
            ``<out>/<sample>/<sample>.trace.txt``, one file per sample.
            **On by default**: it is the record that makes a call auditable
            after the fact, and it costs nothing measurable — 18.5 s against
            19.0 s over five samples, ~160 kB each. A run that did not keep it
            cannot be questioned later without being repeated.

            When ``workers`` is 1 it is *also* streamed to stderr as it
            happens. Parallel workers would interleave the loci of different
            samples into something unreadable, so there the file is the only
            sink.

    Returns:
        One :class:`BatchResult` per entry in input order.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    futures_to_entry: dict[Future[BatchResult], ManifestEntry] = {}
    results_map: dict[str, BatchResult] = {}

    if workers <= 1:
        for entry in entries:
            r = _process_one_sample(
                entry=entry,
                panel=panel,
                out_dir=out_dir,
                reference_fasta=reference_fasta,
                formats=formats,
                min_mapq=min_mapq,
                identity=identity,
                analytical_thresh=analytical_thresh,
                calling_thresh=calling_thresh,
                platform=platform,
                operator=operator,
                run_id=run_id,
                log=log,
                trace=trace,
                trace_live=trace,
            )
            results_map[entry.sample_id] = r
            if progress_callback:
                progress_callback(entry.sample_id)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for entry in entries:
                fut = pool.submit(
                    _process_one_sample,
                    entry=entry,
                    panel=panel,
                    out_dir=out_dir,
                    reference_fasta=reference_fasta,
                    formats=formats,
                    min_mapq=min_mapq,
                    identity=identity,
                    analytical_thresh=analytical_thresh,
                    calling_thresh=calling_thresh,
                    platform=platform,
                    operator=operator,
                    run_id=run_id,
                    log=log,
                    trace=trace,
                )
                futures_to_entry[fut] = entry
            for fut in as_completed(futures_to_entry):
                entry = futures_to_entry[fut]
                try:
                    r = fut.result()
                except Exception as exc:
                    r = BatchResult(
                        sample_id=entry.sample_id,
                        role=entry.role,
                        status="error",
                        error=str(exc),
                    )
                results_map[entry.sample_id] = r
                if progress_callback:
                    progress_callback(entry.sample_id)

    # Return in original manifest order
    results = [results_map[e.sample_id] for e in entries]
    _write_batch_summary(results, panel, out_dir)
    _write_batch_tidy(out_dir, formats)
    _write_cohort_report(results, panel, out_dir, formats)
    return results


@contextmanager
def _trace_sink(
    path: Path | None, *, entry: ManifestEntry, panel: Panel, live: bool = False
) -> Iterator[Callable[[Any], None] | None]:
    """Yield an ``on_trace`` callback that writes one sample's narrative.

    Always to ``path``, one file per sample, so a locus can be re-read
    afterwards. Also to stderr when ``live`` — which the caller sets only for a
    serial run.

    The distinction is interleaving, not volume. With one worker the loci
    arrive in order and watching them is the whole point: the bins, the
    clusters, the aligned sequences and the HP1/HP2 counts are what makes a
    call followable. With several workers the loci of different samples land
    interleaved and the narrative becomes unreadable, so there the file is the
    only sink and ``--log`` is the live channel.

    Yields ``None`` when ``path`` is ``None``, so the caller passes ``on_trace``
    unconditionally and tracing costs nothing when it is off.
    """
    if path is None:
        yield None
        return

    from frontstr.evidence.consensus import poa_backend_name
    from frontstr.interp.naming import default_namer
    from frontstr.panel.stutter_calib import DEFAULT_STUTTER_MODEL
    from frontstr.trace import (
        LocusTrace,
        RunHeader,
        render_header,
        render_locus,
        render_run_summary,
    )
    from frontstr.version import __version__

    namer = default_namer()
    traces: list[LocusTrace] = []
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as fh:
        header = render_header(
            RunHeader(
                inputs=[str(entry.bam)],
                panel_name=panel.name,
                panel_version=panel.version or "",
                n_markers=len(panel.systems),
                consensus_backend=poa_backend_name(),
                naming_markers=sum(1 for s in panel.systems if namer and namer.has_range(s.name)),
                stutter_model=DEFAULT_STUTTER_MODEL.describe(),
                tool_version=__version__,
            )
        )
        fh.write(header + "\n")
        if live:
            print(header, file=sys.stderr, flush=True)

        def emit(locus: LocusTrace) -> None:
            traces.append(locus)
            text = render_locus(locus)
            fh.write(text + "\n\n")
            if live:
                # Flushed per locus: a narrative that appears in blocks when a
                # buffer happens to fill is not something you can follow.
                print(text, end="\n\n", file=sys.stderr, flush=True)

        try:
            yield emit
        finally:
            if traces:
                fh.write(render_run_summary(traces) + "\n")


def _process_one_sample(
    *,
    entry: ManifestEntry,
    panel: Panel,
    out_dir: Path,
    reference_fasta: Path | None,
    formats: frozenset[str],
    min_mapq: int,
    identity: float,
    analytical_thresh: float,
    calling_thresh: float,
    platform: str,
    operator: str | None,
    run_id: str | None,
    log: bool = False,
    trace: bool = False,
    trace_live: bool = False,
) -> BatchResult:
    """Worker function: runs one sample end-to-end and writes output files.

    Must be a module-level function so it is picklable for ProcessPoolExecutor.

    ``log`` and ``trace`` are configured *here*, inside the worker, not in the
    parent: with more than one worker these run in separate processes, and on
    macOS the pool spawns rather than forks, so a logging setup done in the
    parent would not reach them at all.
    """
    import structlog

    from frontstr.exports import (
        write_evidence_csv,
        write_profile_csv,
        write_run_json,
        write_seqs_csv,
    )
    from frontstr.interp import interpret_run
    from frontstr.report import RunContext, build_report, serialize_run

    if log:
        from frontstr.log import configure_logging

        configure_logging(level=logging.DEBUG, console=True)
        # Every line carries its sample. With parallel workers the lines
        # interleave, and an unattributed marker line is worse than no line.
        structlog.contextvars.bind_contextvars(sample=entry.sample_id)

    try:
        sample_dir = out_dir / entry.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        stem = sample_dir / entry.sample_id

        trace_path = stem.with_suffix(".trace.txt") if trace else None
        with _trace_sink(trace_path, entry=entry, panel=panel, live=trace_live) as on_trace:
            marker_results = interpret_run(
                bam=entry.bam,
                panel=panel,
                min_mapq=min_mapq,
                identity_threshold=identity,
                analytical_thresh=analytical_thresh,
                calling_thresh=calling_thresh,
                reference_fasta=reference_fasta,
                on_trace=on_trace,
            )

        context = RunContext(
            sample_name=entry.sample_id,
            panel_name=panel.name,
            panel_version=panel.version,
            bam_path=entry.bam,
            platform=platform,
            operator=operator,
            run_id=run_id,
            reference_build=panel.reference_build,
        )
        payload = serialize_run(marker_results, context)

        written: list[Path] = []
        if "profile" in formats:
            written.append(write_profile_csv(payload, stem.with_suffix(".profile.csv")))
        if "evidence" in formats:
            written.append(write_evidence_csv(payload, stem.with_suffix(".evidence.csv")))
        if "seqs" in formats:
            written.append(write_seqs_csv(payload, stem.with_suffix(".seqs.csv")))
        if "json" in formats:
            written.append(write_run_json(payload, stem.with_suffix(".json"), mode="pretty"))
        if "html" in formats:
            written.append(build_report(marker_results, context, stem.with_suffix(".html")))

        marker_ces = _extract_marker_ces(marker_results)
        return BatchResult(
            sample_id=entry.sample_id,
            role=entry.role,
            status="ok",
            files=written,
            marker_ces=marker_ces,
        )

    except Exception as exc:
        return BatchResult(
            sample_id=entry.sample_id,
            role=entry.role,
            status="error",
            error=f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
        )


def _extract_marker_ces(marker_results: list[MarkerResult]) -> dict[str, str]:
    """Build a {marker_name: allele_string} dict from called alleles.

    The string is e.g. ``"12,14"`` for a het or ``"13"`` for a hom, and ``"X,Y"``
    for AMEL. Undesignatable alleles are rendered as ``"?"``.

    Reads :attr:`Allele.number_label` — the canonical allele number, which comes
    from STRNaming wherever STRNaming has a range for the marker. The raw
    :attr:`Allele.ce` is the length-derived number from before that precedence
    existed; reading it here made ``batch_summary.csv`` the one output still
    reporting pre-STRNaming designations, so a cohort scored on this file
    disagreed with the same run's report, VCF and XLSX at the six markers whose
    ``corr_value`` had been miscalibrated (HG00113 vWA read 13/17, not 14/16).
    """
    out: dict[str, str] = {}
    for r in marker_results:
        parts = [a.number_label or "?" for a in r.alleles_called]
        out[r.marker_name] = ",".join(parts) if parts else ""
    return out


def _write_batch_tidy(out_dir: Path, formats: frozenset[str]) -> list[Path]:
    """Build the cohort tidy dataset from the run JSONs this batch just wrote.

    Reads them back from disk rather than collecting payloads in memory: the
    samples ran in worker processes, and shipping a full payload back through
    pickle for every one of them is a cost paid for no reason when the files
    are already there.

    A no-op when ``json`` was not among the requested formats — without the
    canonical record there is nothing to flatten.
    """
    if "json" not in formats:
        return []
    from frontstr.exports.tidy import load_payloads, write_tidy

    paths = sorted(p for p in out_dir.rglob("*.json") if not p.name.endswith(".min.json"))
    if not paths:
        return []
    try:
        return write_tidy(load_payloads(paths), out_dir)
    except FrontstrError:
        # A malformed run JSON must not sink a batch that otherwise succeeded.
        return []


def _write_cohort_report(
    results: list[BatchResult], panel: Panel, out_dir: Path, formats: frozenset[str]
) -> Path | None:
    """Write ``cohort.html``: every sample at every marker, in one document.

    Only for a real cohort. With one sample the per-sample report already says
    everything, and a second document that repeats it is a file to keep in
    sync for no gain.

    Reads the run JSONs back from disk, the same way the tidy export does, so
    nothing has to travel back from the worker processes. Needs ``json`` for
    the data and ``html`` to be worth linking to, so it is a no-op without
    both.
    """
    if "json" not in formats or "html" not in formats:
        return None
    ok = [r for r in results if r.status == "ok"]
    if len(ok) < 2:
        return None

    from frontstr.exports.tidy import load_payloads
    from frontstr.report.cohort import build_cohort_report

    paths = [out_dir / r.sample_id / f"{r.sample_id}.json" for r in ok]
    existing = [p for p in paths if p.exists()]
    if len(existing) < 2:
        return None

    hrefs = {
        r.sample_id: f"{r.sample_id}/{r.sample_id}.html"
        for r in ok
        if (out_dir / r.sample_id / f"{r.sample_id}.html").exists()
    }
    try:
        return build_cohort_report(
            list(load_payloads(existing)),
            out_dir / "cohort.html",
            report_hrefs=hrefs,
            panel_name=f"{panel.name} {panel.version or ''}".strip(),
        )
    except FrontstrError:
        # A malformed run JSON must not sink a batch that otherwise succeeded,
        # the same rule the tidy export follows.
        return None


def _write_batch_summary(results: list[BatchResult], panel: Panel, out_dir: Path) -> Path:
    """Write ``batch_summary.csv`` to ``out_dir``."""
    marker_names = [s.name for s in panel.systems]
    headers = list(SUMMARY_BASE_HEADERS) + marker_names

    rows: list[dict[str, str]] = []
    for r in results:
        row: dict[str, str] = {
            "sample_id": r.sample_id,
            "role": r.role,
            "status": r.status,
            "error": r.error.splitlines()[0] if r.error else "",
        }
        for m in marker_names:
            row[m] = r.marker_ces.get(m, "") if r.status == "ok" else ""
        rows.append(row)

    summary_path = out_dir / "batch_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return summary_path
