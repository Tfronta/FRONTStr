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
      CTRL_POS/
        ...
      batch_summary.csv
"""

from __future__ import annotations

import csv
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from frontstr.errors import FrontstrError
from frontstr.panel.models import Panel

VALID_ROLES: frozenset[str] = frozenset(
    {"sample", "positive_ctrl", "negative_ctrl", "reagent_blank"}
)

DEFAULT_FORMATS: frozenset[str] = frozenset(
    {"profile", "evidence", "seqs", "json", "html"}
)

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
    status: str                    # "ok" | "error"
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
    if reader.fieldnames is None or "sample_id" not in reader.fieldnames or "bam" not in reader.fieldnames:
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

    Returns:
        One :class:`BatchResult` per entry in input order.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    futures_to_entry: dict = {}
    results_map: dict[str, BatchResult] = {}

    if workers <= 1:
        for entry in entries:
            r = _process_one_sample(
                entry=entry, panel=panel, out_dir=out_dir,
                reference_fasta=reference_fasta, formats=formats,
                min_mapq=min_mapq, identity=identity,
                analytical_thresh=analytical_thresh,
                calling_thresh=calling_thresh,
                platform=platform, operator=operator, run_id=run_id,
            )
            results_map[entry.sample_id] = r
            if progress_callback:
                progress_callback(entry.sample_id)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for entry in entries:
                fut = pool.submit(
                    _process_one_sample,
                    entry=entry, panel=panel, out_dir=out_dir,
                    reference_fasta=reference_fasta, formats=formats,
                    min_mapq=min_mapq, identity=identity,
                    analytical_thresh=analytical_thresh,
                    calling_thresh=calling_thresh,
                    platform=platform, operator=operator, run_id=run_id,
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
    return results


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
) -> BatchResult:
    """Worker function: runs one sample end-to-end and writes output files.

    Must be a module-level function so it is picklable for ProcessPoolExecutor.
    """
    from frontstr.exports import (
        write_evidence_csv,
        write_profile_csv,
        write_run_json,
        write_seqs_csv,
    )
    from frontstr.interp import interpret_run
    from frontstr.report import RunContext, build_report, serialize_run

    try:
        sample_dir = out_dir / entry.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        stem = sample_dir / entry.sample_id

        marker_results = interpret_run(
            bam=entry.bam,
            panel=panel,
            min_mapq=min_mapq,
            identity_threshold=identity,
            analytical_thresh=analytical_thresh,
            calling_thresh=calling_thresh,
            reference_fasta=reference_fasta,
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


def _extract_marker_ces(marker_results: list) -> dict[str, str]:
    """Build a {marker_name: CE_string} dict from called alleles.

    CE_string is e.g. ``"12.0,14.0"`` for het or ``"13.0"`` for hom.
    ``None`` CE values are rendered as ``"?"``.
    """
    out: dict[str, str] = {}
    for r in marker_results:
        parts = []
        for a in r.alleles_called:
            parts.append("?" if a.ce is None else str(a.ce))
        out[r.marker_name] = ",".join(parts) if parts else ""
    return out


def _write_batch_summary(
    results: list[BatchResult], panel: Panel, out_dir: Path
) -> Path:
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
