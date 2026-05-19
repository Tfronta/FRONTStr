"""Tests for :func:`frontstr.report.payload.serialize_run`."""

from __future__ import annotations

import json
from pathlib import Path

from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    MarkerResult,
    TriType,
)
from frontstr.panel.models import System
from frontstr.report.payload import RunContext, serialize_run


def _system(name: str = "TH01") -> System:
    return System(
        name=name, chromosome="chr11", ref_start=2_171_000, ref_end=2_171_050,
        motif="AATG", period=4,
    )


def _allele(idx: int, ce: float, n_reads: int, status: AlleleStatus,
            consensus: str = "AATG") -> Allele:
    return Allele(
        cluster_index=idx, consensus=consensus * int(ce),
        length_bp=len(consensus) * int(ce),
        n_reads_total=n_reads, n_reads_hp1=n_reads // 2, n_reads_hp2=n_reads // 2,
        n_reads_hp_none=0, n_forward=n_reads, n_reverse=0, mean_qual=30.0,
        ce=ce, isfg=f"[{consensus}]{int(ce)}", bp_diff=0, is_deletion=False, status=status,
    )


def _marker_result(name: str, alleles: list[Allele], called: list[Allele],
                   rule: CallRule = CallRule.HETEROZYGOUS,
                   tri: TriType = TriType.NONE,
                   discordant: bool = False) -> MarkerResult:
    return MarkerResult(
        marker_name=name, system=_system(name),
        alleles=alleles, alleles_called=called,
        call_rule=rule, tri_type=tri,
        total_reads=sum(a.n_reads_total for a in alleles),
        discordant=discordant,
    )


def test_serialize_run_minimal() -> None:
    results = [
        _marker_result(
            "TH01",
            [
                _allele(0, 9.0, 60, AlleleStatus.ALLELE),
                _allele(1, 8.0, 55, AlleleStatus.ALLELE),
            ],
            [
                _allele(0, 9.0, 60, AlleleStatus.ALLELE),
                _allele(1, 8.0, 55, AlleleStatus.ALLELE),
            ],
        ),
    ]
    ctx = RunContext(sample_name="S1", panel_name="P1")
    payload = serialize_run(results, ctx)
    assert payload["meta"]["sample_name"] == "S1"
    assert payload["meta"]["app"] == "frontstr"
    assert payload["summary"]["loci_total"] == 1
    assert payload["summary"]["loci_called"] == 1
    assert len(payload["results"]) == 1
    assert len(payload["profile_rows"]) == 1
    row = payload["profile_rows"][0]
    assert row["marker"] == "TH01"
    assert row["allele1_cov"] == 60
    assert row["allele2_cov"] == 55
    assert row["allele3_cov"] is None
    assert row["allele1_ce_label"] == "9"
    assert row["allele1_ce_sort"] == 9
    assert row["allele1_ce_is_kit_ce"] is True
    assert row["status_chip"] == "ok"
    # HP counts must be separate integers, not a combined "hp1/hp2" string.
    assert row["allele1_hp1"] == 30
    assert row["allele1_hp2"] == 30
    assert row["allele2_hp1"] == 27
    assert row["allele2_hp2"] == 27
    assert row["allele3_hp1"] is None
    assert row["allele3_hp2"] is None
    assert "allele1_hp" not in row


def test_profile_row_compound_manual_delta_numeric() -> None:
    """Synthetic compound row: Δ-only allele index from cluster (no YAML ref_ce)."""
    comp = Allele(
        cluster_index=0,
        consensus="A" * 65,
        length_bp=65,
        n_reads_total=33,
        n_reads_hp1=17,
        n_reads_hp2=16,
        n_reads_hp_none=0,
        n_forward=20,
        n_reverse=13,
        mean_qual=30.0,
        ce=None,
        allele_numeric=1.0,
        allele_numeric_source="delta_only",
        isfg="",
        bp_diff=4,
        is_deletion=False,
        status=AlleleStatus.ALLELE,
    )
    results = [
        MarkerResult(
            marker_name="VWA_LIKE",
            system=System(
                name="VWA_LIKE",
                chromosome="chr12",
                ref_start=1,
                ref_end=61,
                motif="TCTA,TCTG",
                period=-1,
                corr_value=0,
            ),
            alleles=[comp],
            alleles_called=[comp],
            call_rule=CallRule.HOMOZYGOUS,
            tri_type=TriType.NONE,
            total_reads=33,
        ),
    ]
    payload = serialize_run(results, RunContext(sample_name="c", panel_name="p"))
    row = payload["profile_rows"][0]
    assert row["allele1_ce"] is None
    assert row["allele1_ce_label"] == "Δ1"
    assert row["allele1_ce_sort"] == 1.0
    assert row["allele1_ce_is_kit_ce"] is False


def test_serialize_run_is_json_safe() -> None:
    results = [
        _marker_result(
            "TPOX",
            [_allele(0, 8.0, 30, AlleleStatus.ALLELE)],
            [_allele(0, 8.0, 30, AlleleStatus.ALLELE)],
            rule=CallRule.HOMOZYGOUS,
        ),
    ]
    payload = serialize_run(results, RunContext(sample_name="X", panel_name="P"))
    # Round-trip via json.dumps to assert there are no non-serializable bits.
    raw = json.dumps(payload, default=str)
    again = json.loads(raw)
    assert again["meta"]["sample_name"] == "X"


def test_serialize_run_tri_and_mixture_flags() -> None:
    cands = [
        _allele(0, 8.0, 30, AlleleStatus.ALLELE),
        _allele(1, 9.0, 28, AlleleStatus.ALLELE),
        _allele(2, 11.0, 25, AlleleStatus.ALLELE),
    ]
    results = [
        _marker_result(
            "TPOX_TRI", cands, cands,
            rule=CallRule.TRIALLELIC_TYPE_II, tri=TriType.TYPE_II_BALANCED,
        ),
        _marker_result(
            "VWA_MIX", cands, cands[:2],
            rule=CallRule.TWO_CALLED_THREE_PRESENT, tri=TriType.MIXTURE_SUSPECTED,
        ),
    ]
    payload = serialize_run(results, RunContext(sample_name="S", panel_name="P"))
    assert payload["summary"]["tri_count"] == 1
    assert payload["summary"]["mixture_count"] == 1
    rows = {r["marker"]: r for r in payload["profile_rows"]}
    assert rows["TPOX_TRI"]["status_chip"] == "tri"
    assert rows["VWA_MIX"]["status_chip"] == "mixture"


def test_serialize_run_dropout_counted() -> None:
    """Markers below the dropout floor must show up in summary.dropouts."""
    small = _marker_result(
        "LOW",
        [_allele(0, 8.0, 10, AlleleStatus.ALLELE)],
        [_allele(0, 8.0, 10, AlleleStatus.ALLELE)],
        rule=CallRule.HOMOZYGOUS,
    )
    no_data = _marker_result("EMPTY", [], [], rule=CallRule.NO_DATA)
    payload = serialize_run([small, no_data], RunContext(sample_name="S", panel_name="P"))
    assert payload["summary"]["dropouts"] == 2  # both are dropouts


def test_serialize_run_computes_bam_hash(tmp_path: Path) -> None:
    bam = tmp_path / "fake.bam"
    bam.write_bytes(b"BAM\x01" * 100)
    ctx = RunContext(sample_name="S", panel_name="P", bam_path=bam)
    payload = serialize_run([], ctx)
    digest = payload["meta"]["bam_sha256"]
    assert isinstance(digest, str)
    assert len(digest) == 64
