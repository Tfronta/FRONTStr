"""CSV exports of a FRONTStr run.

Three flavors, all driven by the same ``payload`` dict that
:func:`frontstr.report.payload.serialize_run` produces:

- :func:`write_profile_csv` — **one row per marker**, wide format with up
  to 3 alleles. Each allele slot carries its own ISFG, CE, integer coverage
  (``alleleN_cov`` plus haplotype split ``alleleN_hp1``/``alleleN_hp2``) and
  consensus sequence (``alleleN_seq``) side by side. Ideal for LIMS import
  and quick spreadsheet review.
- :func:`write_evidence_csv` — **one row per cluster** across all markers,
  long format. Ideal for re-analysis, plotting, and per-allele filtering.
- :func:`write_seqs_csv` — **one row per called allele** with its ISFG
  notation, consensus sequence, length, CE, bp_diff and read count.
  This is the "forensic trail" file: every reported allele's sequence in
  one place.

All files are UTF-8 with ``\\n`` line endings and a stable header order.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from frontstr.interp.isfg import motif_repeat_summary

PROFILE_HEADERS: tuple[str, ...] = (
    "sample",
    "marker",
    "chrom",
    "ref_start",
    "ref_end",
    "motif",
    "call_rule",
    "tri_type",
    "status_chip",
    "total_reads",
    "allele1_isfg",
    "allele1_repeat_summary",
    "allele1_ce",
    "allele1_bp_diff",
    "allele1_cov",
    "allele1_hp1",
    "allele1_hp2",
    "allele1_seq",
    "allele2_isfg",
    "allele2_repeat_summary",
    "allele2_ce",
    "allele2_bp_diff",
    "allele2_cov",
    "allele2_hp1",
    "allele2_hp2",
    "allele2_seq",
    "allele3_isfg",
    "allele3_repeat_summary",
    "allele3_ce",
    "allele3_bp_diff",
    "allele3_cov",
    "allele3_hp1",
    "allele3_hp2",
    "allele3_seq",
)

EVIDENCE_HEADERS: tuple[str, ...] = (
    "sample",
    "marker",
    "cluster_index",
    "status",
    "isfg",
    "isfg_source",
    "motif_repeat_summary",
    "consensus",
    "length_bp",
    "ce",
    "allele_numeric",
    "allele_numeric_source",
    "bp_diff",
    "is_deletion",
    "n_reads_total",
    "n_reads_hp1",
    "n_reads_hp2",
    "n_reads_hp_none",
    "n_forward",
    "n_reverse",
    "mean_qual",
    "expected_stutter",
    "fraction",
)

SEQS_HEADERS: tuple[str, ...] = (
    "sample",
    "marker",
    "allele_index",
    "isfg",
    "isfg_source",
    "motif_repeat_summary",
    "ce",
    "allele_numeric",
    "allele_numeric_source",
    "length_bp",
    "bp_diff",
    "consensus",
    "consensus_method",
    "n_reads_total",
    "n_reads_hp1",
    "n_reads_hp2",
    "status",
)


def write_profile_csv(payload: dict[str, Any], out_path: Path) -> Path:
    """One-row-per-marker forensic profile CSV (wide format)."""
    sample = payload.get("meta", {}).get("sample_name", "")
    rows: list[dict[str, Any]] = []
    for r in payload.get("results", []):
        prow: dict[str, Any] = next(
            (p for p in payload["profile_rows"] if p["marker"] == r["marker_name"]),
            {},
        )
        rows.append(
            {
                "sample": sample,
                "marker": r["marker_name"],
                "chrom": r["system"]["chromosome"],
                "ref_start": r["system"]["ref_start"],
                "ref_end": r["system"]["ref_end"],
                "motif": r["system"]["motif"],
                "call_rule": r["call_rule"],
                "tri_type": r["tri_type"] or "",
                "status_chip": prow.get("status_chip", ""),
                "total_reads": r["total_reads"],
                "allele1_isfg": prow.get("allele1_isfg") or "",
                "allele1_repeat_summary": prow.get("allele1_repeat_summary") or "",
                "allele1_ce": prow.get("allele1_ce_label") or _fmt_number(prow.get("allele1_ce")),
                "allele1_bp_diff": _allele_bp_diff(r, 0),
                "allele1_cov": prow.get("allele1_cov") or "",
                "allele1_hp1": _allele_hp(r, 0, "n_reads_hp1"),
                "allele1_hp2": _allele_hp(r, 0, "n_reads_hp2"),
                "allele1_seq": prow.get("allele1_seq") or "",
                "allele2_isfg": prow.get("allele2_isfg") or "",
                "allele2_repeat_summary": prow.get("allele2_repeat_summary") or "",
                "allele2_ce": prow.get("allele2_ce_label") or _fmt_number(prow.get("allele2_ce")),
                "allele2_bp_diff": _allele_bp_diff(r, 1),
                "allele2_cov": prow.get("allele2_cov") or "",
                "allele2_hp1": _allele_hp(r, 1, "n_reads_hp1"),
                "allele2_hp2": _allele_hp(r, 1, "n_reads_hp2"),
                "allele2_seq": prow.get("allele2_seq") or "",
                "allele3_isfg": prow.get("allele3_isfg") or "",
                "allele3_repeat_summary": prow.get("allele3_repeat_summary") or "",
                "allele3_ce": prow.get("allele3_ce_label") or _fmt_number(prow.get("allele3_ce")),
                "allele3_bp_diff": _allele_bp_diff(r, 2),
                "allele3_cov": prow.get("allele3_cov") or "",
                "allele3_hp1": _allele_hp(r, 2, "n_reads_hp1"),
                "allele3_hp2": _allele_hp(r, 2, "n_reads_hp2"),
                "allele3_seq": prow.get("allele3_seq") or "",
            }
        )
    return _write_csv(out_path, PROFILE_HEADERS, rows)


def write_evidence_csv(payload: dict[str, Any], out_path: Path) -> Path:
    """One-row-per-cluster long-format CSV (all alleles, all statuses)."""
    sample = payload.get("meta", {}).get("sample_name", "")
    rows: list[dict[str, Any]] = []
    for r in payload.get("results", []):
        for a in r.get("alleles", []):
            rows.append(
                {
                    "sample": sample,
                    "marker": r["marker_name"],
                    "cluster_index": a["cluster_index"],
                    "status": a["status"],
                    "isfg": a["isfg"],
                    "isfg_source": a.get("isfg_source", ""),
                    "motif_repeat_summary": motif_repeat_summary(
                        a["consensus"], r["system"]["motif"]
                    ),
                    "consensus": a["consensus"],
                    "length_bp": a["length_bp"],
                    "ce": _fmt_number(a["ce"]),
                    "allele_numeric": _fmt_number(a.get("allele_numeric")),
                    "allele_numeric_source": a.get("allele_numeric_source") or "",
                    "bp_diff": a["bp_diff"],
                    "is_deletion": "true" if a["is_deletion"] else "false",
                    "n_reads_total": a["n_reads_total"],
                    "n_reads_hp1": a["n_reads_hp1"],
                    "n_reads_hp2": a["n_reads_hp2"],
                    "n_reads_hp_none": a["n_reads_hp_none"],
                    "n_forward": a["n_forward"],
                    "n_reverse": a["n_reverse"],
                    "mean_qual": a["mean_qual"],
                    "expected_stutter": a["expected_stutter"],
                    "fraction": a["fraction"],
                }
            )
    return _write_csv(out_path, EVIDENCE_HEADERS, rows)


def write_seqs_csv(payload: dict[str, Any], out_path: Path) -> Path:
    """One-row-per-called-allele forensic sequence trail CSV."""
    sample = payload.get("meta", {}).get("sample_name", "")
    rows: list[dict[str, Any]] = []
    for r in payload.get("results", []):
        for i, a in enumerate(r.get("alleles_called", [])):
            rows.append(
                {
                    "sample": sample,
                    "marker": r["marker_name"],
                    "allele_index": i + 1,
                    "isfg": a["isfg"],
                    "isfg_source": a.get("isfg_source", ""),
                    "motif_repeat_summary": motif_repeat_summary(
                        a["consensus"], r["system"]["motif"]
                    ),
                    "ce": _fmt_number(a["ce"]),
                    "allele_numeric": _fmt_number(a.get("allele_numeric")),
                    "allele_numeric_source": a.get("allele_numeric_source") or "",
                    "length_bp": a["length_bp"],
                    "bp_diff": a["bp_diff"],
                    "consensus": a["consensus"],
                    "consensus_method": a.get("consensus_method", ""),
                    "n_reads_total": a["n_reads_total"],
                    "n_reads_hp1": a["n_reads_hp1"],
                    "n_reads_hp2": a["n_reads_hp2"],
                    "status": a["status"],
                }
            )
    return _write_csv(out_path, SEQS_HEADERS, rows)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _allele_bp_diff(result: dict[str, Any], i: int) -> Any:
    called = result.get("alleles_called", [])
    return called[i]["bp_diff"] if i < len(called) else ""


def _allele_hp(result: dict[str, Any], i: int, key: str) -> Any:
    called = result.get("alleles_called", [])
    return called[i][key] if i < len(called) else ""


def _fmt_number(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        # Trim trailing zeros for prettier CSVs, but keep at least one decimal.
        s = f"{value:.4f}".rstrip("0")
        return s.rstrip(".") if s.endswith(".") else s
    return str(value)


def _write_csv(out_path: Path, headers: tuple[str, ...], rows: list[dict[str, Any]]) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return out_path
