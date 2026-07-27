"""Tests for :mod:`frontstr.audit` — process log, audit record, integrity hash."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from frontstr.audit import (
    AUDIT_SCHEMA_VERSION,
    AuditRecord,
    InputFile,
    build_audit_record,
    configure_logging,
    file_sha256,
    get_logger,
)
from frontstr.interp.models import (
    CallRule,
    Flag,
    FlagCode,
    FlagSeverity,
    MarkerResult,
    TriType,
)
from frontstr.interp.qc import QcThresholds
from frontstr.panel.models import System


@pytest.fixture(autouse=True)
def _restore_quiet_logging():
    """Leave logging as importing the library found it."""
    yield
    configure_logging(None, level=logging.WARNING)


def _result(name: str = "TH01", flags: list[Flag] | None = None) -> MarkerResult:
    return MarkerResult(
        marker_name=name,
        system=System(
            name=name, chromosome="chr11", ref_start=1, ref_end=400, motif="AATG", period=4
        ),
        alleles=[],
        alleles_called=[],
        call_rule=CallRule.HETEROZYGOUS,
        tri_type=TriType.NONE,
        total_reads=40,
        flags=flags or [],
    )


# ---------------------------------------------------------------------------
# Process log
# ---------------------------------------------------------------------------


def test_process_log_is_one_json_object_per_line(tmp_path: Path) -> None:
    log_path = tmp_path / "run.jsonl"
    configure_logging(log_path, level=logging.INFO)
    log = get_logger("test")
    log.info("first.event", marker="TH01", reads=40)
    log.info("second.event", ok=True)
    logging.shutdown()

    lines = [json.loads(x) for x in log_path.read_text().splitlines() if x.strip()]
    assert [x["event"] for x in lines] == ["first.event", "second.event"]
    assert lines[0]["marker"] == "TH01"
    assert lines[0]["level"] == "info"
    assert "timestamp" in lines[0]


def test_debug_events_are_dropped_at_info_level(tmp_path: Path) -> None:
    log_path = tmp_path / "run.jsonl"
    configure_logging(log_path, level=logging.INFO)
    get_logger("test").debug("noisy.event")
    logging.shutdown()
    assert log_path.read_text().strip() == ""


def test_library_import_does_not_configure_logging(tmp_path: Path) -> None:
    """A library that writes to someone else's stdout is a broken library."""
    configure_logging(None, level=logging.WARNING)
    log_path = tmp_path / "run.jsonl"
    handler = logging.FileHandler(log_path, encoding="utf-8")
    logging.basicConfig(handlers=[handler], level=logging.DEBUG, force=True)
    get_logger("test").info("should.not.appear")
    logging.shutdown()
    assert "should.not.appear" not in log_path.read_text()


# ---------------------------------------------------------------------------
# Audit record
# ---------------------------------------------------------------------------


def test_record_counts_flags_by_code_and_severity() -> None:
    results = [
        _result("TH01", [Flag.of(FlagCode.LOW_COVERAGE, "x")]),
        _result("vWA", [Flag.of(FlagCode.LOW_COVERAGE, "y"), Flag.of(FlagCode.ISOALLELE, "z")]),
    ]
    record = build_audit_record(results)
    assert record.flag_counts == {"isoallele": 1, "low_coverage": 2}
    assert record.severity_counts == {"info": 1, "warn": 2}


def test_markers_needing_review_excludes_info_only_markers() -> None:
    """An INFO flag is context, not a reason to pull the marker for review."""
    results = [
        _result("TH01", [Flag.of(FlagCode.ISOALLELE, "info only")]),
        _result("vWA", [Flag.of(FlagCode.LOW_COVERAGE, "a warning")]),
    ]
    assert build_audit_record(results).markers_needing_review == ["vWA"]


def test_flags_checked_lists_every_known_code() -> None:
    """The guarantee that makes a clean report meaningful.

    Without it, a code missing from the counts is ambiguous between "checked
    and not found" and "never looked at".
    """
    record = build_audit_record([_result()])
    assert set(record.flags_checked) == {c.value for c in FlagCode}
    assert record.flag_counts == {}


def test_record_captures_the_resolved_environment() -> None:
    """Not what the caller intended — what the run actually used."""
    record = build_audit_record([_result()])
    assert record.poa_backend, "must state a backend, even if 'none'"
    assert record.stutter_model_version
    assert record.stutter_model_protocol == "wgs_pcr_free"
    assert record.schema_version == AUDIT_SCHEMA_VERSION


def test_thresholds_are_recorded() -> None:
    record = build_audit_record(
        [_result()],
        qc_thresholds=QcThresholds(low_coverage_reads=17),
        analytical_thresh=0.03,
        calling_thresh=0.15,
    )
    assert record.qc_thresholds.low_coverage_reads == 17
    assert (record.analytical_thresh, record.calling_thresh) == (0.03, 0.15)


def test_allele_level_flags_are_counted_too() -> None:
    from frontstr.interp.models import Allele, AlleleStatus

    allele = Allele(
        cluster_index=0,
        consensus="A",
        length_bp=1,
        n_reads_total=5,
        n_reads_hp1=0,
        n_reads_hp2=0,
        n_reads_hp_none=5,
        n_forward=3,
        n_reverse=2,
        mean_qual=30.0,
        ce=None,
        isfg="",
        bp_diff=0,
        is_deletion=False,
        status=AlleleStatus.ALLELE,
        flags=[Flag.of(FlagCode.STRAND_BIAS, "on the allele")],
    )
    result = _result()
    result.alleles = [allele]
    record = build_audit_record([result])
    assert record.flag_counts == {"strand_bias": 1}
    assert record.markers_needing_review == ["TH01"]


# ---------------------------------------------------------------------------
# Integrity
# ---------------------------------------------------------------------------


def test_built_record_is_sealed_and_verifies() -> None:
    record = build_audit_record([_result()])
    assert record.integrity_sha256 is not None
    assert record.verify()


def test_editing_the_record_breaks_the_seal() -> None:
    record = build_audit_record([_result()])
    tampered = record.model_copy(update={"poa_backend": "something else"})
    assert not tampered.verify()


def test_unsealed_record_does_not_verify() -> None:
    """No hash means no claim; it must not read as 'verified'."""
    assert not AuditRecord().verify()


def test_seal_is_stable_across_equal_records() -> None:
    """Field insertion order must not change the hash."""
    a = AuditRecord(poa_backend="poa_spoa", flag_counts={"a": 1, "b": 2}).sealed()
    b = AuditRecord(flag_counts={"b": 2, "a": 1}, poa_backend="poa_spoa").sealed()
    assert a.integrity_sha256 == b.integrity_sha256


def test_seal_covers_nested_input_hashes() -> None:
    base = AuditRecord(inputs=[InputFile(role="bam", path="a.bam", sha256="deadbeef")])
    sealed = base.sealed()
    swapped = sealed.model_copy(
        update={"inputs": [InputFile(role="bam", path="a.bam", sha256="cafebabe")]}
    )
    assert not swapped.verify()


def test_file_sha256_matches_hashlib(tmp_path: Path) -> None:
    import hashlib

    path = tmp_path / "x.bin"
    payload = b"forensic" * 1000
    path.write_bytes(payload)
    assert file_sha256(path) == hashlib.sha256(payload).hexdigest()


def test_severity_enum_values_used_by_review_triage() -> None:
    """Guards the string values the triage logic branches on."""
    assert {s.value for s in FlagSeverity} >= {"info", "warn", "error"}
