"""Cohort-scale tidy dataset.

Every other export describes one sample. This one describes a cohort: N runs
collapsed into a single long table, one row per **sample × marker × allele**,
which is the shape an analysis actually wants. Grouping up to the genotype
level is a one-line ``group_by``; going the other way is not possible, so the
finest grain is the one worth storing.

It is built from the canonical run JSONs rather than from live results. That
matters for three reasons: batch runs its samples in worker processes and the
payloads never come back; a dataset can be rebuilt at any time without
re-running the caller; and runs from different batches — different days,
different panel versions — can be combined into one table.

Dropouts are rows too
---------------------

A marker with no called allele emits one row with ``allele_index = 0``,
``called = false`` and null allele fields. Omitting it would make a dropout
indistinguishable from a marker that was never in the panel — which is exactly
the distinction a concordance study needs.

Why the run configuration rides on every row
--------------------------------------------

``panel_version``, ``poa_backend`` and ``stutter_model`` are repeated per row.
That is redundant in the CSV and nearly free in Parquet, which
dictionary-encodes them. It is there because a 150-sample benchmark is not
collected in one afternoon: when half the cohort was called under one stutter
model and half under another, the dataset should say so rather than leave you
to reconstruct it from timestamps.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from frontstr.errors import FrontstrError

#: Column order. Identity first, then the allele, then coverage, then the
#: marker and run context. ``consensus`` is last because it is long and would
#: otherwise push everything useful off the right of a spreadsheet.
TIDY_COLUMNS: tuple[str, ...] = (
    "sample",
    "marker",
    "allele_index",
    "called",
    "allele_number",
    "allele_label",
    "number_method",
    "number_is_absolute",
    "isfg",
    "iso_suffix",
    "iso_match",
    "length_bp",
    "reads",
    "reads_hp1",
    "reads_hp2",
    "reads_absorbed",
    "allele_fraction",
    "consensus_method",
    "call_rule",
    "locus_reads",
    "n_alleles_called",
    "marker_flags",
    "needs_review",
    "panel_name",
    "panel_version",
    "poa_backend",
    "stutter_model",
    "run_id",
    "consensus",
)

#: Explicit Parquet types. Inferring from Python rows misreads a column that is
#: all-null in the first chunk, and a benchmark that silently types
#: ``allele_number`` as string will compare "10" < "9".
_PARQUET_TYPES: dict[str, str] = {
    "allele_index": "int32",
    "called": "bool",
    "allele_number": "float64",
    "number_is_absolute": "bool",
    "length_bp": "int32",
    "reads": "int32",
    "reads_hp1": "int32",
    "reads_hp2": "int32",
    "reads_absorbed": "int32",
    "allele_fraction": "float64",
    "locus_reads": "int32",
    "n_alleles_called": "int32",
    "needs_review": "bool",
}


def load_payloads(paths: Iterable[Path]) -> Iterator[dict[str, Any]]:
    """Read canonical run JSONs, yielding one payload each.

    Raises:
        FrontstrError: If a file is not readable JSON or is not a FRONTStr run.
    """
    for path in paths:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise FrontstrError(f"Cannot read run JSON {path}: {exc}") from exc
        if not isinstance(payload, dict) or "results" not in payload:
            raise FrontstrError(
                f"{path} is not a FRONTStr run JSON (no 'results' key). "
                "Pass the .json written by `frontstr export --formats json`."
            )
        yield payload


def _run_context(payload: dict[str, Any]) -> dict[str, Any]:
    meta = payload.get("meta", {})
    audit = payload.get("audit", {})
    return {
        "sample": meta.get("sample_name"),
        "panel_name": meta.get("panel_name"),
        "panel_version": meta.get("panel_version"),
        "poa_backend": audit.get("poa_backend"),
        "stutter_model": audit.get("stutter_model_version"),
        "run_id": meta.get("run_id"),
    }


def _marker_context(marker: dict[str, Any]) -> dict[str, Any]:
    flags = marker.get("flags", [])
    return {
        "marker": marker["marker_name"],
        "call_rule": marker["call_rule"],
        "locus_reads": marker["total_reads"],
        "n_alleles_called": len(marker.get("alleles_called", [])),
        "marker_flags": ",".join(f["code"] for f in flags) or None,
        "needs_review": any(f["severity"] in ("warn", "error") for f in flags),
    }


def _allele_fields(allele: dict[str, Any], index: int) -> dict[str, Any]:
    iso = allele.get("iso") or {}
    suffix = iso.get("suffix")
    return {
        "allele_index": index,
        "called": True,
        "allele_number": allele.get("number"),
        "allele_label": allele.get("number_label"),
        "number_method": allele.get("number_method"),
        "number_is_absolute": allele.get("number_is_absolute"),
        "isfg": allele.get("isfg") or None,
        "iso_suffix": suffix,
        "iso_match": iso.get("match_type") if suffix else None,
        "length_bp": allele.get("length_bp"),
        "reads": allele.get("n_reads_total"),
        "reads_hp1": allele.get("n_reads_hp1"),
        "reads_hp2": allele.get("n_reads_hp2"),
        "reads_absorbed": allele.get("n_reads_absorbed", 0),
        "allele_fraction": allele.get("fraction"),
        "consensus_method": allele.get("consensus_method"),
        "consensus": allele.get("consensus") or None,
    }


_EMPTY_ALLELE: dict[str, Any] = {
    "allele_index": 0,
    "called": False,
    "allele_number": None,
    "allele_label": None,
    "number_method": None,
    "number_is_absolute": None,
    "isfg": None,
    "iso_suffix": None,
    "iso_match": None,
    "length_bp": None,
    "reads": None,
    "reads_hp1": None,
    "reads_hp2": None,
    "reads_absorbed": None,
    "allele_fraction": None,
    "consensus_method": None,
    "consensus": None,
}


def build_tidy_rows(payloads: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Flatten run payloads into one row per sample × marker × allele.

    Markers with no called allele contribute a single ``called = False`` row so
    a dropout stays visible in the dataset.
    """
    rows: list[dict[str, Any]] = []
    for payload in payloads:
        run = _run_context(payload)
        for marker in payload.get("results", []):
            context = {**run, **_marker_context(marker)}
            called = marker.get("alleles_called", [])
            if not called:
                rows.append({**context, **_EMPTY_ALLELE})
                continue
            for i, allele in enumerate(called, start=1):
                rows.append({**context, **_allele_fields(allele, i)})
    return rows


def write_tidy_csv(rows: list[dict[str, Any]], out_path: Path) -> Path:
    """Write the tidy dataset as CSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIDY_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in TIDY_COLUMNS})
    return out_path


def parquet_available() -> bool:
    """True when pyarrow is importable."""
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        return False
    return True


def write_tidy_parquet(rows: list[dict[str, Any]], out_path: Path) -> Path:
    """Write the tidy dataset as Parquet, with an explicit schema.

    Raises:
        FrontstrError: If pyarrow is not installed.
    """
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise FrontstrError(
            "Parquet export needs pyarrow: pip install 'frontstr[parquet]'"
        ) from exc

    type_map = {
        "int32": pa.int32(),
        "float64": pa.float64(),
        "bool": pa.bool_(),
    }
    schema = pa.schema(
        [
            pa.field(name, type_map.get(_PARQUET_TYPES.get(name, ""), pa.string()))
            for name in TIDY_COLUMNS
        ]
    )
    columns = {name: [row.get(name) for row in rows] for name in TIDY_COLUMNS}
    table = pa.Table.from_pydict(columns, schema=schema)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path, compression="zstd")
    return out_path


def write_tidy(
    payloads: Iterable[dict[str, Any]],
    out_dir: Path,
    *,
    stem: str = "cohort_tidy",
) -> list[Path]:
    """Build and write the tidy dataset as both CSV and Parquet.

    Parquet is skipped — not failed — when pyarrow is absent, so a cohort
    export still produces something usable on a machine without it. The caller
    is expected to surface that; see the CLI's warning.

    Returns:
        The paths written, CSV first.
    """
    rows = build_tidy_rows(payloads)
    written = [write_tidy_csv(rows, out_dir / f"{stem}.csv")]
    if parquet_available():
        written.append(write_tidy_parquet(rows, out_dir / f"{stem}.parquet"))
    return written
