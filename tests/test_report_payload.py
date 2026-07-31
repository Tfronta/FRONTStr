"""Tests for :func:`frontstr.report.payload.serialize_run`."""

from __future__ import annotations

import json
from pathlib import Path

from frontstr.interp.flags import derive_marker_flags
from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    Flag,
    FlagCode,
    IsoAllele,
    MarkerResult,
    TriType,
)
from frontstr.interp.qc import QcThresholds
from frontstr.panel.models import System
from frontstr.report.payload import RunContext, serialize_run


def _system(name: str = "TH01") -> System:
    return System(
        name=name,
        chromosome="chr11",
        ref_start=2_171_000,
        ref_end=2_171_050,
        motif="AATG",
        period=4,
    )


def _allele(
    idx: int, ce: float, n_reads: int, status: AlleleStatus, consensus: str = "AATG"
) -> Allele:
    return Allele(
        cluster_index=idx,
        consensus=consensus * int(ce),
        length_bp=len(consensus) * int(ce),
        n_reads_total=n_reads,
        n_reads_hp1=n_reads // 2,
        n_reads_hp2=n_reads // 2,
        n_reads_hp_none=0,
        n_forward=n_reads,
        n_reverse=0,
        mean_qual=30.0,
        ce=ce,
        isfg=f"[{consensus}]{int(ce)}",
        bp_diff=0,
        is_deletion=False,
        status=status,
    )


def _marker_result(
    name: str,
    alleles: list[Allele],
    called: list[Allele],
    rule: CallRule = CallRule.HETEROZYGOUS,
    tri: TriType = TriType.NONE,
    discordant: bool = False,
) -> MarkerResult:
    return MarkerResult(
        marker_name=name,
        system=_system(name),
        alleles=alleles,
        alleles_called=called,
        call_rule=rule,
        tri_type=tri,
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
    # CE table headline: combined genotype + no isoalleles here.
    assert row["genotype"] == "9 / 8"
    assert row["has_iso"] is False
    # Sequencing-based flat table: one row per called allele.
    assert len(payload["seq_rows"]) == 2
    s0 = payload["seq_rows"][0]
    assert s0["marker"] == "TH01"
    assert s0["number"] == "9"
    assert s0["iso"] == ""
    assert s0["isfg"] == "[AATG]9"
    assert s0["consensus"] == "AATG" * 9
    assert s0["n_reads_total"] == 60
    # HP counts must be separate integers, not a combined "hp1/hp2" string.
    assert row["allele1_hp1"] == 30
    assert row["allele1_hp2"] == 30
    assert row["allele2_hp1"] == 27
    assert row["allele2_hp2"] == 27
    assert row["allele3_hp1"] is None
    assert row["allele3_hp2"] is None
    assert "allele1_hp" not in row


def test_iso_alleles_flag_and_designation() -> None:
    """Catalog suffix → iso designation in seq_rows + has_iso on the CE row."""
    a1 = Allele(
        cluster_index=0,
        consensus="TCTA" * 15,
        length_bp=60,
        n_reads_total=19,
        n_reads_hp1=18,
        n_reads_hp2=1,
        n_reads_hp_none=0,
        n_forward=12,
        n_reverse=7,
        mean_qual=30.0,
        ce=15.0,
        isfg="[TCTA]15",
        bp_diff=0,
        is_deletion=False,
        allele_numeric=15.0,
        allele_numeric_source="period_ce",
        status=AlleleStatus.ALLELE,
        iso=IsoAllele(suffix="a", match_type="exact", distance=0, source="STRSeq"),
    )
    a2 = Allele(
        cluster_index=1,
        consensus="TCTG" + "TCTA" * 14,
        length_bp=60,
        n_reads_total=17,
        n_reads_hp1=1,
        n_reads_hp2=16,
        n_reads_hp_none=0,
        n_forward=9,
        n_reverse=8,
        mean_qual=30.0,
        ce=15.0,
        isfg="TCTG [TCTA]14",
        bp_diff=0,
        is_deletion=False,
        allele_numeric=15.0,
        allele_numeric_source="period_ce",
        status=AlleleStatus.ALLELE,
        iso=IsoAllele(suffix="b", match_type="exact", distance=0, source="STRSeq"),
    )
    result = _marker_result("D3S1358", [a1, a2], [a1, a2])
    derive_marker_flags(result)
    payload = serialize_run([result], RunContext(sample_name="c", panel_name="p"))
    row = payload["profile_rows"][0]
    assert row["genotype"] == "15 / 15"
    assert row["has_iso"] is True
    iso_designations = {s["iso"] for s in payload["seq_rows"]}
    assert iso_designations == {"15a", "15b"}


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


def test_profile_row_compound_shows_bracket_count_not_delta() -> None:
    """Compound marker with a bracket-derived ``ce`` shows the absolute repeat
    count as the allele number, never the relative Δ offset."""
    comp = Allele(
        cluster_index=0,
        consensus="TCTA" * 13,
        length_bp=261,
        n_reads_total=16,
        n_reads_hp1=8,
        n_reads_hp2=8,
        n_reads_hp_none=0,
        n_forward=10,
        n_reverse=6,
        mean_qual=30.0,
        ce=13.0,
        allele_numeric=-3.0,
        allele_numeric_source="delta_only",
        isfg="[TCTA]13",
        bp_diff=-12,
        is_deletion=False,
        status=AlleleStatus.ALLELE,
    )
    results = [
        MarkerResult(
            marker_name="vWA",
            system=System(
                name="vWA",
                chromosome="chr12",
                ref_start=5_983_877,
                ref_end=5_984_149,
                motif="TCTA,TCTG",
                period=-1,
                corr_value=8,
            ),
            alleles=[comp],
            alleles_called=[comp],
            call_rule=CallRule.HOMOZYGOUS,
            tri_type=TriType.NONE,
            total_reads=16,
        ),
    ]
    payload = serialize_run(results, RunContext(sample_name="c", panel_name="p"))
    row = payload["profile_rows"][0]
    # Absolute repeat count from the sequence, not the Δ-3 offset.
    assert row["allele1_ce_label"] == "13"
    assert row["allele1_ce_sort"] == 13.0
    assert row["allele1_ce_is_kit_ce"] is True
    assert "Δ" not in row["allele1_ce_label"]


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
            "TPOX_TRI",
            cands,
            cands,
            rule=CallRule.TRIALLELIC_TYPE_II,
            tri=TriType.TYPE_II_BALANCED,
        ),
        _marker_result(
            "VWA_MIX",
            cands,
            cands[:2],
            rule=CallRule.TWO_CALLED_THREE_PRESENT,
            tri=TriType.MIXTURE_SUSPECTED,
        ),
    ]
    payload = serialize_run(results, RunContext(sample_name="S", panel_name="P"))
    assert payload["summary"]["tri_count"] == 1
    assert payload["summary"]["mixture_count"] == 1
    rows = {r["marker"]: r for r in payload["profile_rows"]}
    assert rows["TPOX_TRI"]["status_chip"] == "tri"
    assert rows["VWA_MIX"]["status_chip"] == "mixture"


def test_low_coverage_kpi_counts_flags_not_a_threshold_of_its_own() -> None:
    """The cover page reports the flags the caller raised, and nothing else.

    It used to re-derive the condition from a report-local floor of 30 against
    the window's spanning depth, and label the result "drop-outs". On the
    bundled sample that read 12 next to "25 called", none of the 12 having
    dropped out. A marker is thin here only if it carries LOW_COVERAGE.
    """
    thin = _marker_result(
        "THIN",
        [_allele(0, 8.0, 10, AlleleStatus.ALLELE)],
        [_allele(0, 8.0, 10, AlleleStatus.ALLELE)],
        rule=CallRule.HOMOZYGOUS,
    )
    thin.flags.append(Flag.of(FlagCode.LOW_COVERAGE, "called on 10 reads"))

    # Under the retired floor of 30 and unflagged: must not be counted.
    unflagged = _marker_result(
        "SHALLOW_BUT_FINE",
        [_allele(0, 8.0, 25, AlleleStatus.ALLELE)],
        [_allele(0, 8.0, 25, AlleleStatus.ALLELE)],
        rule=CallRule.HOMOZYGOUS,
    )

    payload = serialize_run([thin, unflagged], RunContext(sample_name="S", panel_name="P"))
    assert payload["summary"]["low_coverage"] == 1
    assert "dropouts" not in payload["summary"]


def test_low_coverage_floor_comes_from_the_qc_policy() -> None:
    """The floor in the payload is the one the calls were made under."""
    ctx = RunContext(
        sample_name="S", panel_name="P", qc_thresholds=QcThresholds(low_coverage_reads=12)
    )
    payload = serialize_run([], ctx)
    assert payload["meta"]["low_coverage_floor"] == 12
    assert "dropout_floor" not in payload["meta"]


def test_serialize_run_computes_bam_hash(tmp_path: Path) -> None:
    bam = tmp_path / "fake.bam"
    bam.write_bytes(b"BAM\x01" * 100)
    ctx = RunContext(sample_name="S", panel_name="P", bam_path=bam)
    payload = serialize_run([], ctx)
    digest = payload["meta"]["bam_sha256"]
    assert isinstance(digest, str)
    assert len(digest) == 64
