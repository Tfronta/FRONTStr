"""Single-source serializer of an interpreted run.

Every consumer of FRONTStr's results — HTML report, CSV/XLSX exports, the
future REST API — calls :func:`serialize_run` so the data is identical
everywhere. The output is a plain JSON-ready ``dict`` (no Pydantic, no
dataclasses) so it can also be inlined into the HTML report as
``<script type="application/json">``.

"""

from __future__ import annotations

import hashlib
from collections import Counter
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from frontstr.audit import InputFile, build_audit_record
from frontstr.interp.isfg import motif_repeat_summary
from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    FlagCode,
    MarkerResult,
    TriType,
)
from frontstr.interp.qc import QcThresholds
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
    reference_build: str = "GRCh38"
    platform: str = "ont"
    operator: str | None = None
    run_id: str | None = None
    started_at: datetime = field(default_factory=lambda: datetime.now(UTC))
    #: The exact command line, for the report's provenance section. Populated
    #: by the CLI from ``sys.argv``; empty for library callers.
    pipeline_argv: list[str] = field(default_factory=list)
    #: Every parameter actually in force, defaults included — ``pipeline_argv``
    #: alone is misleading, because the values that decide a call are usually
    #: the ones nobody typed. Rows from
    #: :meth:`frontstr.params.RunParameters.as_audit_rows`, so each carries its
    #: default, whether it was changed, and where the default came from.
    effective_params: list[dict[str, Any]] = field(default_factory=list)
    #: The panel's extraction windows as BED lines. See
    #: :func:`frontstr.panel.bed.panel_bed_lines`.
    panel_bed: list[str] = field(default_factory=list)
    dropout_floor: int = DEFAULT_DROPOUT_FLOOR
    #: QC policy applied during interpretation, carried through so the audit
    #: record states the thresholds the calls were actually made under.
    qc_thresholds: QcThresholds = field(default_factory=QcThresholds)


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
            "pipeline_argv": context.pipeline_argv,
            "effective_params": context.effective_params,
            "panel_bed": context.panel_bed,
            "dropout_floor": context.dropout_floor,
        },
        "summary": summary,
        "qc": qc,
        "audit": _build_audit(results, context),
        "profile_rows": [_profile_row(r, context.sample_name) for r in results],
        "seq_rows": _seq_rows(results),
        "results": serialized_results,
    }
    payload["strhub"] = build_strhub_projection(payload)
    return payload


def _build_audit(results: list[MarkerResult], context: RunContext) -> dict[str, Any]:
    """Assemble the audit record from the run's inputs and its results.

    Placed in the canonical payload rather than only in the HTML: the audit
    trail has to travel with the data, or it is not an audit trail.
    """
    inputs = [
        InputFile(role=role, path=str(path), sha256=sha)
        for role, path, sha in (
            ("bam", context.bam_path, context.bam_sha256),
            ("panel", None, context.panel_sha256),
        )
        if path is not None or sha is not None
    ]
    first = results[0] if results else None
    record = build_audit_record(
        results,
        inputs=inputs,
        qc_thresholds=context.qc_thresholds,
        analytical_thresh=first.analytical_thresh if first else None,
        calling_thresh=first.calling_thresh if first else None,
    )
    return record.model_dump(mode="json")


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
            "strand": r.system.strand,
            "marker_type": r.system.marker_type,
            "allow_triallelic": r.system.allow_triallelic,
        },
        "call_rule": r.call_rule.value,
        "tri_type": r.tri_type.value,
        "total_reads": r.total_reads,
        "analytical_thresh": r.analytical_thresh,
        "calling_thresh": r.calling_thresh,
        "flags": [f.model_dump(mode="json") for f in r.flags],
        "alleles": [
            _serialize_allele(a, r.total_reads, r.system.motif, r.system.strand) for a in r.alleles
        ],
        "alleles_called": [
            _serialize_allele(a, r.total_reads, r.system.motif, r.system.strand)
            for a in r.alleles_called
        ],
    }
    marker_dict["ngs_panel"] = build_ngs_panel(marker_dict)
    return marker_dict


def _serialize_allele(a: Allele, total_reads: int, motif: str, strand: str = "+") -> dict[str, Any]:
    return {
        "cluster_index": a.cluster_index,
        "consensus": a.consensus,
        "length_bp": a.length_bp,
        "ce": a.ce,
        "number": a.number,
        "number_method": a.number_method,
        "number_label": a.number_label,
        "number_is_absolute": a.number_is_absolute,
        "allele_numeric": a.allele_numeric,
        "allele_numeric_source": a.allele_numeric_source,
        # The canonical bracketed string: STRNaming's when it has a range for
        # this marker, the legacy full-window scan otherwise. One string per
        # allele across every view — see Allele.repeat_label.
        "isfg": a.repeat_label,
        "isfg_source": a.repeat_label_source,
        # The raw window scan, kept so nothing is lost: it spans the whole
        # extraction window rather than the standard reporting range.
        "isfg_window": a.isfg,
        "motif_repeat_summary": motif_repeat_summary(a.consensus, motif, strand=strand),
        "bp_diff": a.bp_diff,
        "is_deletion": a.is_deletion,
        "consensus_method": a.consensus_method,
        "n_reads_total": a.n_reads_total,
        "n_reads_hp1": a.n_reads_hp1,
        "n_reads_hp2": a.n_reads_hp2,
        "n_reads_hp_none": a.n_reads_hp_none,
        "n_forward": a.n_forward,
        "n_reverse": a.n_reverse,
        "mean_qual": round(a.mean_qual, 2),
        "n_reads_absorbed": a.n_reads_absorbed,
        "expected_stutter": round(a.expected_stutter, 3),
        "status": a.status.value,
        "fraction": round(a.fraction(total_reads), 4),
        "iso": a.iso.model_dump(mode="json"),
        "flags": [f.model_dump(mode="json") for f in a.flags],
    }


def _trim_ce_display(ce: float) -> str:
    """Pretty CE string from a forensic CE value."""
    x = round(float(ce), 4)
    if abs(x - int(x)) < 1e-9:
        return str(round(x))
    return f"{x:.10f}".rstrip("0").rstrip(".")


def _format_allele_number(a: Allele) -> tuple[float | None, str | None, bool]:
    """Table cell for an allele: ``(sort key, label, is_absolute_number)``.

    A thin adapter over the model. The number, its label and whether it is a
    real absolute allele number are all decided on :class:`Allele` so that no
    view can render an allele differently from another.
    """
    return (a.number, a.number_label or None, a.number_is_absolute)


def _profile_row(r: MarkerResult, sample_name: str = "") -> dict[str, Any]:
    """Wide row for the profile table: 1 marker, up to 3 alleles.

    Carries the **flags** as well as the coarse status chip. They used to live
    only inside the expandable per-locus cards at the bottom of the report,
    which meant a reviewer scanning the profile table could not see that a
    locus was flagged at all — the XLSX export had been doing this correctly
    (tinted rows plus a QC sheet) while the HTML did not.
    """
    called = list(r.alleles_called)[:3]
    row: dict[str, Any] = {
        "sample": sample_name,
        "marker": r.marker_name,
        "call_rule": r.call_rule.value,
        "tri_type": r.tri_type.value,
        "total_reads": r.total_reads,
        # The two numbers the CLI already reports separately. `total_reads` is
        # the denominator every fraction threshold is measured against, so it
        # stays; but it is the wrong figure to *show* as the locus coverage —
        # see MarkerResult.called_reads.
        "called_reads": r.called_reads,
        "discarded_reads": r.discarded_reads,
        "status_chip": _status_chip(r),
        "allele_balance": r.allele_balance,
        # `short` for the chip, `code` and `message` for its tooltip and for
        # the legend under the table. All three from FlagCode, so the chips
        # and the legend cannot describe different things.
        "flags": [
            {
                "code": f.code.value,
                "short": f.code.short,
                "severity": f.severity.value,
                "message": f.message,
            }
            for f in r.flags
        ],
        "worst_severity": _worst_severity(r),
    }
    for i in range(3):
        slot = called[i] if i < len(called) else None
        if slot is not None:
            row[f"allele{i + 1}_isfg"] = slot.repeat_label
            row[f"allele{i + 1}_repeat_summary"] = motif_repeat_summary(
                slot.consensus, r.system.motif, strand=r.system.strand
            )
            row[f"allele{i + 1}_ce"] = slot.ce
            ce_sort, ce_label, ce_is_kit = _format_allele_number(slot)
            row[f"allele{i + 1}_ce_sort"] = ce_sort
            row[f"allele{i + 1}_ce_label"] = ce_label
            row[f"allele{i + 1}_ce_is_kit_ce"] = ce_is_kit
            row[f"allele{i + 1}_cov"] = slot.n_reads_total
            row[f"allele{i + 1}_hp1"] = slot.n_reads_hp1
            row[f"allele{i + 1}_hp2"] = slot.n_reads_hp2
            row[f"allele{i + 1}_seq"] = slot.consensus
        else:
            row[f"allele{i + 1}_isfg"] = None
            row[f"allele{i + 1}_repeat_summary"] = None
            row[f"allele{i + 1}_ce"] = None
            row[f"allele{i + 1}_ce_sort"] = None
            row[f"allele{i + 1}_ce_label"] = None
            row[f"allele{i + 1}_ce_is_kit_ce"] = None
            row[f"allele{i + 1}_cov"] = None
            row[f"allele{i + 1}_hp1"] = None
            row[f"allele{i + 1}_hp2"] = None
            row[f"allele{i + 1}_seq"] = None

    # Genotype string for the CE table headline (numeric, e.g. "12 / 14").
    labels = [row[f"allele{i + 1}_ce_label"] for i in range(3) if row[f"allele{i + 1}_ce_label"]]
    row["genotype"] = " / ".join(labels) if labels else "\u2013"

    # Isoallele presence is decided in the model layer (interp.flags); the CE
    # table only reads the resulting marker flag.
    row["has_iso"] = any(f.code == FlagCode.ISOALLELE for f in r.flags)
    return row


def _worst_severity(r: MarkerResult) -> str:
    """``error`` | ``warn`` | ``info`` | ``""`` — drives the row tint."""
    severities = {f.severity.value for f in r.flags}
    for level in ("error", "warn", "info"):
        if level in severities:
            return level
    return ""


def _allele_number_label(a: Allele) -> str | None:
    """The allele-number cell for one allele (identical to the CE table's)."""
    return a.number_label or None


def _seq_rows(results: list[MarkerResult]) -> list[dict[str, Any]]:
    """Flat sequencing-based table: one row per called allele, all markers.

    This is the NGS differential view (ISFG / iso-allele / full sequence),
    kept separate from the CE-based genotype table. The raw consensus travels
    here so consumers can offer copy / FASTA without bloating the CE table.
    """
    rows: list[dict[str, Any]] = []
    for r in results:
        for i, a in enumerate(r.alleles_called):
            number = _allele_number_label(a)
            iso = f"{number}{a.iso.suffix}" if (number and a.iso.suffix) else ""
            rows.append(
                {
                    "marker": r.marker_name,
                    "allele_index": i + 1,
                    "number": number or "\u2013",
                    "iso": iso,
                    "n_reads_total": a.n_reads_total,
                    "n_reads_hp1": a.n_reads_hp1,
                    "n_reads_hp2": a.n_reads_hp2,
                    "isfg": a.repeat_label,
                    "isfg_source": a.repeat_label_source,
                    "motif_repeat_summary": motif_repeat_summary(
                        a.consensus, r.system.motif, strand=r.system.strand
                    ),
                    "length_bp": a.length_bp,
                    "consensus": a.consensus,
                    "status": a.status.value,
                }
            )
    return rows


def _status_chip(r: MarkerResult) -> str:
    """Top-level status badge: ok | low | tri | mixture | no_data."""
    if r.call_rule == CallRule.NO_DATA:
        return "no_data"
    if r.tri_type == TriType.MIXTURE_SUSPECTED:
        return "mixture"
    if r.tri_type in (TriType.TYPE_I_UNBALANCED, TriType.TYPE_II_BALANCED):
        return "tri"
    return "ok"


def _compute_summary(results: list[MarkerResult], dropout_floor: int) -> dict[str, Any]:
    """KPI numbers shown on the cover page."""
    loci_total = len(results)
    loci_called = sum(1 for r in results if r.call_rule != CallRule.NO_DATA)
    tri_count = sum(
        1 for r in results if r.tri_type in (TriType.TYPE_I_UNBALANCED, TriType.TYPE_II_BALANCED)
    )
    mixture_count = sum(1 for r in results if r.tri_type == TriType.MIXTURE_SUSPECTED)
    dropouts = sum(
        1 for r in results if 0 < r.total_reads < dropout_floor or r.call_rule == CallRule.NO_DATA
    )
    return {
        "loci_total": loci_total,
        "loci_called": loci_called,
        "tri_count": tri_count,
        "mixture_count": mixture_count,
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
        "status_breakdown": [{"status": s, "count": c} for s, c in sorted(statuses.items())],
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
