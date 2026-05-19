"""Single-source serializer of an interpreted run.

Every consumer of FRONTStr's results — HTML report, CSV/XLSX exports, the
future REST API — calls :func:`serialize_run` so the data is identical
everywhere. The output is a plain JSON-ready ``dict`` (no Pydantic, no
dataclasses) so it can also be inlined into the HTML report as
``<script type="application/json">``.

The shape follows plan-longtr-improved.md §17.4.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from frontstr.interp.isfg import motif_repeat_summary
from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    MarkerResult,
    TriType,
)
from frontstr.report.ngs_display import build_ngs_panel, build_strhub_projection
from frontstr.version import __version__

DEFAULT_DROPOUT_FLOOR = 30


@dataclass(slots=True)
class RunContext:
    """Non-result metadata about the run.

    Hash fields are optional; if ``bam_path`` is provided the BAM hash is
    computed automatically. Everything else flows through verbatim.
    """

    sample_name: str
    panel_name: str
    panel_version: str = ""
    panel_sha256: str | None = None
    bam_path: Path | None = None
    bam_sha256: str | None = None
    longtr_vcf_path: Path | None = None
    longtr_vcf_sha256: str | None = None
    longtr_version: str | None = None
    reference_build: str = "GRCh38"
    platform: str = "ont"
    operator: str | None = None
    run_id: str | None = None
    started_at: datetime = field(default_factory=lambda: datetime.now(UTC))
    pipeline_argv: list[str] = field(default_factory=list)
    dropout_floor: int = DEFAULT_DROPOUT_FLOOR


def serialize_run(
    results: list[MarkerResult],
    context: RunContext,
) -> dict[str, Any]:
    """Turn the interpretation output + context into a JSON-ready dict.

    Args:
        results: One :class:`MarkerResult` per marker, in panel order.
        context: Run-level metadata (sample, panel, hashes, operator).

    Returns:
        Dict with top-level keys ``meta`` / ``summary`` / ``qc`` /
        ``profile_rows`` / ``results`` / ``strhub``. Each marker in
        ``results`` includes ``ngs_panel`` for HTML / tooling. Safe to
        ``json.dumps`` directly.
    """
    if context.bam_path is not None and context.bam_sha256 is None:
        context.bam_sha256 = _file_sha256(context.bam_path)
    if context.longtr_vcf_path is not None and context.longtr_vcf_sha256 is None:
        context.longtr_vcf_sha256 = _file_sha256(context.longtr_vcf_path)

    serialized_results = [_serialize_marker(r) for r in results]
    summary = _compute_summary(results, context.dropout_floor)
    qc = _compute_qc(results)

    payload: dict[str, Any] = {
        "meta": {
            "app": "frontstr",
            "app_version": __version__,
            "sample_name": context.sample_name,
            "operator": context.operator,
            "run_id": context.run_id,
            "started_at": context.started_at.isoformat(),
            "platform": context.platform,
            "reference_build": context.reference_build,
            "panel_name": context.panel_name,
            "panel_version": context.panel_version,
            "panel_sha256": context.panel_sha256,
            "bam_path": str(context.bam_path) if context.bam_path else None,
            "bam_sha256": context.bam_sha256,
            "longtr_vcf_path": (
                str(context.longtr_vcf_path) if context.longtr_vcf_path else None
            ),
            "longtr_vcf_sha256": context.longtr_vcf_sha256,
            "longtr_version": context.longtr_version,
            "pipeline_argv": context.pipeline_argv,
            "dropout_floor": context.dropout_floor,
        },
        "summary": summary,
        "qc": qc,
        "profile_rows": [_profile_row(r) for r in results],
        "results": serialized_results,
    }
    payload["strhub"] = build_strhub_projection(payload)
    return payload


def _serialize_marker(r: MarkerResult) -> dict[str, Any]:
    marker_dict: dict[str, Any] = {
        "marker_name": r.marker_name,
        "system": {
            "name": r.system.name,
            "codis_name": r.system.codis_name,
            "chromosome": r.system.chromosome,
            "ref_start": r.system.ref_start,
            "ref_end": r.system.ref_end,
            "motif": r.system.motif,
            "period": r.system.period,
            "corr_value": r.system.corr_value,
            "reference_ce": r.system.reference_ce,
            "allele_bp_step": r.system.allele_bp_step,
            "category": r.system.category,
            "allow_triallelic": r.system.allow_triallelic,
        },
        "call_rule": r.call_rule.value,
        "tri_type": r.tri_type.value,
        "total_reads": r.total_reads,
        "analytical_thresh": r.analytical_thresh,
        "calling_thresh": r.calling_thresh,
        "discordant": r.discordant,
        "warnings": list(r.warnings),
        "alleles": [_serialize_allele(a, r.total_reads, r.system.motif) for a in r.alleles],
        "alleles_called": [
            _serialize_allele(a, r.total_reads, r.system.motif) for a in r.alleles_called
        ],
        "longtr": _serialize_longtr(r),
    }
    marker_dict["ngs_panel"] = build_ngs_panel(marker_dict)
    return marker_dict


def _serialize_allele(a: Allele, total_reads: int, motif: str) -> dict[str, Any]:
    return {
        "cluster_index": a.cluster_index,
        "consensus": a.consensus,
        "length_bp": a.length_bp,
        "ce": a.ce,
        "allele_numeric": a.allele_numeric,
        "allele_numeric_source": a.allele_numeric_source,
        "isfg": a.isfg,
        "motif_repeat_summary": motif_repeat_summary(a.consensus, motif),
        "bp_diff": a.bp_diff,
        "is_deletion": a.is_deletion,
        "n_reads_total": a.n_reads_total,
        "n_reads_hp1": a.n_reads_hp1,
        "n_reads_hp2": a.n_reads_hp2,
        "n_reads_hp_none": a.n_reads_hp_none,
        "n_forward": a.n_forward,
        "n_reverse": a.n_reverse,
        "mean_qual": round(a.mean_qual, 2),
        "expected_stutter": round(a.expected_stutter, 3),
        "status": a.status.value,
        "longtr_match": a.longtr_match,
        "longtr_inexact": a.longtr_inexact,
        "longtr_bp_diff": a.longtr_bp_diff,
        "fraction": round(a.fraction(total_reads), 4),
    }


def _serialize_longtr(r: MarkerResult) -> dict[str, Any] | None:
    if r.longtr_result is None:
        return None
    lt = r.longtr_result
    sample_call = next(iter(lt.samples.values()), None)
    return {
        "marker_name": lt.marker_name,
        "chrom": lt.chrom,
        "pos": lt.pos,
        "motif": lt.motif,
        "period": lt.period,
        "alleles": [
            {
                "sequence": a.sequence,
                "bp_diff": a.bp_diff,
                "inexact": a.inexact,
                "is_deletion": a.is_deletion,
            }
            for a in lt.alleles
        ],
        "gt_indices": list(sample_call.gt_indices) if sample_call and sample_call.gt_indices else None,
        "posterior": (
            round(sample_call.posterior, 4)
            if sample_call and sample_call.posterior is not None
            else None
        ),
        "depth": sample_call.depth if sample_call else 0,
        "pdp_hp1": sample_call.pdp_hp1 if sample_call else 0,
        "pdp_hp2": sample_call.pdp_hp2 if sample_call else 0,
    }


def _trim_ce_display(ce: float) -> str:
    """Pretty CE string from a forensic CE value."""
    x = round(float(ce), 4)
    if abs(x - int(x)) < 1e-9:
        return str(int(round(x)))
    return f"{x:.10f}".rstrip("0").rstrip(".")


def _profile_allele_ce_sort_and_label(
    ce: float | None, length_bp: int
) -> tuple[float | None, str | None, bool]:
    """Legacy CE / TR-bp cell when :attr:`Allele.allele_numeric` is absent."""
    if ce is not None:
        x = round(float(ce), 4)
        return (x, _trim_ce_display(ce), True)
    if length_bp > 0:
        return (float(length_bp), f"{length_bp} bp TR", False)
    return (None, None, False)


def _profile_forensic_display(
    *,
    allele_numeric: float | None,
    allele_src: str,
    ce: float | None,
    length_bp: int,
) -> tuple[float | None, str | None, bool]:
    """(sort key, table cell, confident_absolute_scale).

    Prefer :class:`Allele.allele_numeric`; fall back to CE / bp TR strings.
    """
    if (
        allele_numeric is not None
        and allele_src
        and allele_src not in {"unavailable", "deletion"}
    ):
        base = _trim_ce_display(float(allele_numeric))
        if allele_src == "delta_only" and abs(float(allele_numeric)) > 1e-9:
            label = f"Δ{base}"
        else:
            label = base
        confident = allele_src in {"period_ce", "reference_offset"}
        return float(allele_numeric), label, confident
    return _profile_allele_ce_sort_and_label(ce, length_bp)


def _profile_row(r: MarkerResult) -> dict[str, Any]:
    """Wide row for the profile table: 1 marker, up to 3 alleles."""
    called = list(r.alleles_called)[:3]
    row: dict[str, Any] = {
        "marker": r.marker_name,
        "call_rule": r.call_rule.value,
        "tri_type": r.tri_type.value,
        "total_reads": r.total_reads,
        "discordant": r.discordant,
        "status_chip": _status_chip(r),
    }
    for i in range(3):
        slot = called[i] if i < len(called) else None
        if slot is not None:
            row[f"allele{i + 1}_isfg"] = slot.isfg
            row[f"allele{i + 1}_repeat_summary"] = motif_repeat_summary(
                slot.consensus, r.system.motif
            )
            row[f"allele{i + 1}_ce"] = slot.ce
            ce_sort, ce_label, ce_is_kit = _profile_forensic_display(
                allele_numeric=slot.allele_numeric,
                allele_src=slot.allele_numeric_source,
                ce=slot.ce,
                length_bp=slot.length_bp,
            )
            row[f"allele{i + 1}_ce_sort"] = ce_sort
            row[f"allele{i + 1}_ce_label"] = ce_label
            row[f"allele{i + 1}_ce_is_kit_ce"] = ce_is_kit
            row[f"allele{i + 1}_cov"] = slot.n_reads_total
            row[f"allele{i + 1}_hp"] = f"{slot.n_reads_hp1}/{slot.n_reads_hp2}"
        else:
            row[f"allele{i + 1}_isfg"] = None
            row[f"allele{i + 1}_repeat_summary"] = None
            row[f"allele{i + 1}_ce"] = None
            row[f"allele{i + 1}_ce_sort"] = None
            row[f"allele{i + 1}_ce_label"] = None
            row[f"allele{i + 1}_ce_is_kit_ce"] = None
            row[f"allele{i + 1}_cov"] = None
            row[f"allele{i + 1}_hp"] = None
    return row


def _status_chip(r: MarkerResult) -> str:
    """Top-level status badge: ok | low | tri | mixture | discordant | no_data."""
    if r.call_rule == CallRule.NO_DATA:
        return "no_data"
    if r.tri_type == TriType.MIXTURE_SUSPECTED:
        return "mixture"
    if r.tri_type in (TriType.TYPE_I_UNBALANCED, TriType.TYPE_II_BALANCED):
        return "tri"
    if r.discordant:
        return "discordant"
    return "ok"


def _compute_summary(results: list[MarkerResult], dropout_floor: int) -> dict[str, Any]:
    """KPI numbers shown on the cover page."""
    loci_total = len(results)
    loci_called = sum(1 for r in results if r.call_rule != CallRule.NO_DATA)
    tri_count = sum(
        1 for r in results if r.tri_type in (TriType.TYPE_I_UNBALANCED, TriType.TYPE_II_BALANCED)
    )
    mixture_count = sum(1 for r in results if r.tri_type == TriType.MIXTURE_SUSPECTED)
    discordant = sum(1 for r in results if r.discordant)
    dropouts = sum(
        1
        for r in results
        if 0 < r.total_reads < dropout_floor
        or r.call_rule == CallRule.NO_DATA
    )
    return {
        "loci_total": loci_total,
        "loci_called": loci_called,
        "tri_count": tri_count,
        "mixture_count": mixture_count,
        "discordant_count": discordant,
        "dropouts": dropouts,
    }


def _compute_qc(results: list[MarkerResult]) -> dict[str, Any]:
    """Aggregated QC metrics for the dashboard."""
    covs = [r.total_reads for r in results]
    statuses: Counter[str] = Counter()
    for r in results:
        for a in r.alleles:
            statuses[a.status.value] += 1
    coverage_table = [
        {"marker": r.marker_name, "coverage": r.total_reads, "chip": _status_chip(r)}
        for r in results
    ]
    return {
        "mean_coverage": (sum(covs) / len(covs)) if covs else 0.0,
        "min_coverage": min(covs) if covs else 0,
        "max_coverage": max(covs) if covs else 0,
        "coverage_table": coverage_table,
        "status_breakdown": [
            {"status": s, "count": c} for s, c in sorted(statuses.items())
        ],
    }


def _file_sha256(path: Path, chunk: int = 1 << 20) -> str:
    """Stream a file through SHA-256 (chunk size 1 MiB)."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for buf in iter(lambda: fh.read(chunk), b""):
            h.update(buf)
    return h.hexdigest()


# Re-export for convenience
__all__ = [
    "DEFAULT_DROPOUT_FLOOR",
    "AlleleStatus",
    "RunContext",
    "serialize_run",
]
