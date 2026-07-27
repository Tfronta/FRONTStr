"""Tests for CSV / JSON exports.

We build the same payload used by the report tests and assert that the CSVs
have stable headers, the right number of rows, correct values for tricky
cases (deletions, tri-allelic, missing CE), and that JSON round-trips.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from frontstr.errors import FrontstrError
from frontstr.exports import (
    EVIDENCE_HEADERS,
    PROFILE_HEADERS,
    SEQS_HEADERS,
    write_evidence_csv,
    write_profile_csv,
    write_run_json,
    write_seqs_csv,
)
from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    MarkerResult,
    TriType,
)
from frontstr.panel.models import System
from frontstr.report.payload import RunContext, serialize_run


def _allele(
    idx: int,
    ce: float | None,
    cov: int,
    status: AlleleStatus,
    *,
    consensus: str = "AATGAATGAATG",
    is_del: bool = False,
) -> Allele:
    return Allele(
        cluster_index=idx,
        consensus="" if is_del else consensus,
        length_bp=0 if is_del else len(consensus),
        n_reads_total=cov,
        n_reads_hp1=cov // 2,
        n_reads_hp2=cov - cov // 2,
        n_reads_hp_none=0,
        n_forward=cov,
        n_reverse=0,
        mean_qual=30.0,
        ce=ce,
        isfg=f"[AATG]{int(ce)}" if ce is not None else "",
        bp_diff=-12 if is_del else 0,
        is_deletion=is_del,
        status=status,
    )


def _system(name: str, tri: bool = False) -> System:
    return System(
        name=name,
        chromosome="chr2",
        ref_start=1_489_651,
        ref_end=1_489_684,
        motif="AATG",
        period=4,
        allow_triallelic=tri,
        tri_balanced_thr=0.5 if tri else None,
    )


@pytest.fixture
def small_payload() -> dict:
    th01 = MarkerResult(
        marker_name="TH01",
        system=_system("TH01"),
        alleles=[
            _allele(0, 9.0, 60, AlleleStatus.ALLELE),
            _allele(1, 8.0, 55, AlleleStatus.ALLELE),
            _allele(2, 7.0, 4, AlleleStatus.STUTTER),
        ],
        alleles_called=[
            _allele(0, 9.0, 60, AlleleStatus.ALLELE),
            _allele(1, 8.0, 55, AlleleStatus.ALLELE),
        ],
        call_rule=CallRule.HETEROZYGOUS,
        tri_type=TriType.NONE,
        total_reads=119,
    )
    tpox = MarkerResult(
        marker_name="TPOX",
        system=_system("TPOX", tri=True),
        alleles=[
            _allele(0, 8.0, 30, AlleleStatus.ALLELE),
            _allele(1, 9.0, 28, AlleleStatus.ALLELE),
            _allele(2, 11.0, 26, AlleleStatus.ALLELE),
        ],
        alleles_called=[
            _allele(0, 8.0, 30, AlleleStatus.ALLELE),
            _allele(1, 9.0, 28, AlleleStatus.ALLELE),
            _allele(2, 11.0, 26, AlleleStatus.ALLELE),
        ],
        call_rule=CallRule.TRIALLELIC_TYPE_II,
        tri_type=TriType.TYPE_II_BALANCED,
        total_reads=84,
    )
    del_marker = MarkerResult(
        marker_name="DEL_LOCUS",
        system=_system("DEL_LOCUS"),
        alleles=[_allele(0, None, 40, AlleleStatus.DELETION, is_del=True)],
        alleles_called=[_allele(0, None, 40, AlleleStatus.DELETION, is_del=True)],
        call_rule=CallRule.HOMOZYGOUS,
        tri_type=TriType.NONE,
        total_reads=40,
    )
    return serialize_run(
        [th01, tpox, del_marker],
        RunContext(sample_name="S001", panel_name="demo", panel_version="0.1"),
    )


def test_profile_csv_headers_stable(small_payload: dict, tmp_path: Path) -> None:
    out = write_profile_csv(small_payload, tmp_path / "p.csv")
    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert rows[0]["sample"] == "S001"
    assert tuple(rows[0].keys()) == PROFILE_HEADERS
    by_marker = {r["marker"]: r for r in rows}
    th01 = by_marker["TH01"]
    assert th01["allele1_cov"] == "60"
    assert th01["allele2_cov"] == "55"
    assert th01["allele3_cov"] == ""
    # Per-allele consensus sequence must travel alongside coverage.
    assert th01["allele1_seq"] != ""
    assert th01["allele2_seq"] != ""
    assert th01["allele3_seq"] == ""
    assert th01["call_rule"] == "heterozygous"
    assert th01["tri_type"] == ""
    tpox = by_marker["TPOX"]
    assert tpox["allele3_cov"] == "26"
    assert tpox["tri_type"] == "tri_II_balanced"
    assert tpox["status_chip"] == "tri"


def test_profile_csv_handles_deletion(small_payload: dict, tmp_path: Path) -> None:
    out = write_profile_csv(small_payload, tmp_path / "p.csv")
    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    del_row = next(r for r in rows if r["marker"] == "DEL_LOCUS")
    # Deletion has no CE and no ISFG, but coverage and bp_diff must be present
    assert del_row["allele1_ce"] == ""
    assert del_row["allele1_isfg"] == ""
    assert del_row["allele1_cov"] == "40"
    assert del_row["allele1_bp_diff"] == "-12"


def test_evidence_csv_one_row_per_cluster(small_payload: dict, tmp_path: Path) -> None:
    out = write_evidence_csv(small_payload, tmp_path / "e.csv")
    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert tuple(rows[0].keys()) == EVIDENCE_HEADERS
    # 3 + 3 + 1 = 7 clusters
    assert len(rows) == 7
    statuses = {r["status"] for r in rows}
    assert {"allele", "stutter", "deletion"}.issubset(statuses)
    # Stutter row must have fraction > 0 but n_reads small
    stutter = next(r for r in rows if r["status"] == "stutter")
    assert int(stutter["n_reads_total"]) == 4
    deletion = next(r for r in rows if r["status"] == "deletion")
    assert deletion["is_deletion"] == "true"


def test_seqs_csv_only_called_alleles(small_payload: dict, tmp_path: Path) -> None:
    out = write_seqs_csv(small_payload, tmp_path / "s.csv")
    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert tuple(rows[0].keys()) == SEQS_HEADERS
    # 2 + 3 + 1 = 6 called alleles across the 3 markers
    assert len(rows) == 6
    # Each TPOX row must have allele_index 1..3
    tpox_indices = sorted(int(r["allele_index"]) for r in rows if r["marker"] == "TPOX")
    assert tpox_indices == [1, 2, 3]
    # Consensus strings present for the non-deletion alleles
    th01_first = next(r for r in rows if r["marker"] == "TH01" and r["allele_index"] == "1")
    assert th01_first["consensus"] != ""


def test_json_export_round_trip(small_payload: dict, tmp_path: Path) -> None:
    pretty = write_run_json(small_payload, tmp_path / "pretty.json", mode="pretty")
    compact = write_run_json(small_payload, tmp_path / "compact.json", mode="compact")
    # Pretty must have newlines + indent; compact must not
    assert "\n  " in pretty.read_text(encoding="utf-8")
    assert "\n " not in compact.read_text(encoding="utf-8").rstrip("\n")
    # Both parse back to the same dict
    a = json.loads(pretty.read_text(encoding="utf-8"))
    b = json.loads(compact.read_text(encoding="utf-8"))
    assert a == b
    assert a["meta"]["sample_name"] == "S001"
    assert len(a["results"]) == 3
    # Compact is meaningfully smaller
    assert compact.stat().st_size < pretty.stat().st_size


def test_json_export_unknown_mode(small_payload: dict, tmp_path: Path) -> None:
    with pytest.raises(FrontstrError, match="Unknown JSON mode"):
        write_run_json(small_payload, tmp_path / "x.json", mode="binary")  # type: ignore[arg-type]


def test_csv_empty_payload(tmp_path: Path) -> None:
    """An empty payload must still write a valid CSV (header only)."""
    empty = serialize_run([], RunContext(sample_name="E", panel_name="P"))
    p = write_profile_csv(empty, tmp_path / "empty.csv")
    lines = p.read_text(encoding="utf-8").splitlines()
    assert lines[0] == ",".join(PROFILE_HEADERS)
    assert len(lines) == 1
