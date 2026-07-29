"""Tests for :func:`frontstr.report.html.build_report`.

We parse the generated HTML with lxml in lenient mode to ensure it is at
least well-formed enough for a browser, and assert that the key data points
made it into the document.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import lxml.html

from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    MarkerResult,
    TriType,
)
from frontstr.interp.qc import derive_run_qc_flags
from frontstr.panel.models import System
from frontstr.report import RunContext, build_report


def _allele(
    idx: int, ce: float, cov: int, status: AlleleStatus, consensus_motif: str = "AATG"
) -> Allele:
    return Allele(
        cluster_index=idx,
        consensus=consensus_motif * int(ce),
        length_bp=len(consensus_motif) * int(ce),
        n_reads_total=cov,
        n_reads_hp1=cov // 2,
        n_reads_hp2=cov - cov // 2,
        n_reads_hp_none=0,
        n_forward=cov,
        n_reverse=0,
        mean_qual=30.0,
        ce=ce,
        isfg=f"[{consensus_motif}]{int(ce)}",
        bp_diff=0,
        is_deletion=False,
        status=status,
    )


def _make_results() -> list[MarkerResult]:
    th01_system = System(
        name="TH01",
        chromosome="chr11",
        ref_start=2_171_000,
        ref_end=2_171_050,
        motif="AATG",
        period=4,
    )
    tpox_system = System(
        name="TPOX",
        chromosome="chr2",
        ref_start=1_489_651,
        ref_end=1_489_684,
        motif="AATG",
        period=4,
        allow_triallelic=True,
        tri_balanced_thr=0.5,
    )
    th01_alleles = [
        _allele(0, 9.0, 60, AlleleStatus.ALLELE),
        _allele(1, 8.0, 55, AlleleStatus.ALLELE),
        _allele(2, 7.0, 6, AlleleStatus.STUTTER),
    ]
    tpox_alleles = [
        _allele(0, 8.0, 30, AlleleStatus.ALLELE),
        _allele(1, 9.0, 28, AlleleStatus.ALLELE),
        _allele(2, 11.0, 26, AlleleStatus.ALLELE),
    ]
    return [
        MarkerResult(
            marker_name="TH01",
            system=th01_system,
            alleles=th01_alleles,
            alleles_called=th01_alleles[:2],
            call_rule=CallRule.HETEROZYGOUS,
            tri_type=TriType.NONE,
            total_reads=sum(a.n_reads_total for a in th01_alleles),
        ),
        MarkerResult(
            marker_name="TPOX",
            system=tpox_system,
            alleles=tpox_alleles,
            alleles_called=tpox_alleles,
            call_rule=CallRule.TRIALLELIC_TYPE_II,
            tri_type=TriType.TYPE_II_BALANCED,
            total_reads=sum(a.n_reads_total for a in tpox_alleles),
        ),
    ]


def test_build_report_writes_valid_html(tmp_path: Path) -> None:
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(
            sample_name="S001", panel_name="forensic-panel", panel_version="0.1", operator="J. Diaz"
        ),
        out,
    )
    assert out.exists()
    html = out.read_text(encoding="utf-8")
    assert html.lstrip().lower().startswith("<!doctype html>")
    doc = lxml.html.fromstring(html)

    # Sample name in <title> and <h1>
    title = doc.find(".//title").text
    assert "S001" in title
    h1 = doc.find(".//h1")
    assert h1 is not None and "S001" in (h1.text or "")

    # Required sections present (incl. the sequencing-based table)
    for sid in ("cover", "profile", "sequences", "qc", "loci", "raw"):
        assert doc.get_element_by_id(sid) is not None

    # CE table = genotype headline; sequencing table = one row per called allele
    assert ">Genotype<" in html
    seqtable = doc.find(".//table[@class='profile seqtable']")
    assert seqtable is not None
    seq_body_rows = seqtable.findall(".//tbody/tr")
    assert len(seq_body_rows) == 5  # TH01 het (2) + TPOX tri (3)
    # Full consensus reaches the sequencing table behind a copy control
    assert "data-copy=" in html
    # CE table marker links over to the sequencing table, not the locus detail
    assert 'class="marker-link"' in html

    # Both markers rendered as <details>
    th01 = doc.get_element_by_id("locus-TH01")
    tpox = doc.get_element_by_id("locus-TPOX")
    assert th01.tag == "details"
    assert tpox.tag == "details"
    # Tri-allelic locus must be open by default
    assert tpox.get("open") is not None


def test_build_report_embeds_json_payload(tmp_path: Path) -> None:
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(sample_name="S2", panel_name="P", operator=None),
        out,
    )
    html = out.read_text(encoding="utf-8")
    doc = lxml.html.fromstring(html)
    blob = doc.get_element_by_id("run-data")
    assert blob is not None
    data = json.loads(blob.text)
    assert data["meta"]["sample_name"] == "S2"
    assert {r["marker_name"] for r in data["results"]} == {"TH01", "TPOX"}
    assert data["strhub"]["schema"] == "strhub.ngs_panel/v1"
    assert len(data["strhub"]["markers"]) == 2


def test_build_report_stamps_self_hash(tmp_path: Path) -> None:
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(sample_name="S3", panel_name="P"),
        out,
    )
    html = out.read_text(encoding="utf-8")
    matches = re.findall(r"[0-9a-f]{64}", html)
    assert matches, "expected at least one SHA-256 hash in the report footer"
    # Hash should appear both inline and in footer (two slots)
    sha_in_footer = re.search(r'id="report-hash-footer">([0-9a-f]{64})<', html)
    assert sha_in_footer is not None


def test_build_report_includes_ngs_panel(tmp_path: Path) -> None:
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(sample_name="S-NGS", panel_name="P"),
        out,
    )
    html = out.read_text(encoding="utf-8")
    assert "Next-Generation Sequencing Analysis" in html
    assert 'class="ngs-panel"' in html
    assert "(reads)" in html
    assert "Forensic audit" in html
    assert "sequence-scroll" in html


def test_profile_tri_columns_only_when_tri_chip(tmp_path: Path) -> None:
    """Allele 3 columns appear iff some marker has status_chip tri."""
    out_tri = tmp_path / "tri.html"
    build_report(
        _make_results(),
        RunContext(sample_name="T", panel_name="P"),
        out_tri,
    )
    html_tri = out_tri.read_text(encoding="utf-8")
    assert ">Allele 3</th>" in html_tri

    th01_only = [
        MarkerResult(
            marker_name="TH01",
            system=System(
                name="TH01",
                chromosome="chr11",
                ref_start=2_171_000,
                ref_end=2_171_050,
                motif="AATG",
                period=4,
            ),
            alleles=[
                _allele(0, 9.0, 60, AlleleStatus.ALLELE),
                _allele(1, 8.0, 55, AlleleStatus.ALLELE),
            ],
            alleles_called=[
                _allele(0, 9.0, 60, AlleleStatus.ALLELE),
                _allele(1, 8.0, 55, AlleleStatus.ALLELE),
            ],
            call_rule=CallRule.HETEROZYGOUS,
            tri_type=TriType.NONE,
            total_reads=115,
        ),
    ]
    out_no = tmp_path / "no_tri.html"
    build_report(
        th01_only,
        RunContext(sample_name="N", panel_name="P"),
        out_no,
    )
    html_no = out_no.read_text(encoding="utf-8")
    assert ">Allele 3</th>" not in html_no


def test_build_report_is_deterministic(tmp_path: Path) -> None:
    """Same inputs must produce byte-identical reports (sha-stable)."""
    results = _make_results()
    ctx = RunContext(sample_name="DET", panel_name="P", panel_version="1")
    a = tmp_path / "a.html"
    b = tmp_path / "b.html"
    # Freeze the started_at so the meta block matches across runs
    ctx.started_at = ctx.started_at.replace(microsecond=0)
    build_report(results, ctx, a)
    build_report(results, ctx, b)
    assert hashlib.sha256(a.read_bytes()).digest() == hashlib.sha256(b.read_bytes()).digest()


def test_build_report_includes_chips_for_special_calls(tmp_path: Path) -> None:
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(sample_name="S", panel_name="P"),
        out,
    )
    html = out.read_text(encoding="utf-8")
    assert "tri_II_balanced" in html
    assert "TPOX" in html
    # Status chips for both alleles
    assert html.count('class="chip ok"') > 0
    assert "triallelic_type_II" in html


# ---------------------------------------------------------------------------
# Profile table: everything on one row (plan item 2 — flags were invisible here)
# ---------------------------------------------------------------------------


def _flagged_results() -> list[MarkerResult]:
    """TH01 heterozygous but unevenly covered, so QC has something to show."""
    alleles = [
        _allele(0, 9.0, 20, AlleleStatus.ALLELE),
        _allele(1, 8.0, 6, AlleleStatus.ALLELE),
    ]
    r = MarkerResult(
        marker_name="TH01",
        system=System(
            name="TH01",
            chromosome="chr11",
            ref_start=2_171_000,
            ref_end=2_171_050,
            motif="AATG",
            period=4,
        ),
        alleles=alleles,
        alleles_called=alleles,
        call_rule=CallRule.HETEROZYGOUS,
        tri_type=TriType.NONE,
        total_reads=26,
    )
    derive_run_qc_flags([r])
    return [r]


def _profile_table(html: str):
    return lxml.html.fromstring(html).find(".//table[@class='profile profile-wide']")


def test_profile_row_carries_sample_alleles_coverage_and_sequences(tmp_path: Path) -> None:
    """One row per locus must answer the whole question without a second view."""
    out = tmp_path / "r.html"
    build_report(_make_results(), RunContext(sample_name="S-ROW", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")

    headers = [th.text_content().strip() for th in _profile_table(html).findall(".//thead/tr/th")]
    joined = " | ".join(headers)
    for wanted in ("Sample", "Marker", "Allele 1", "Reads 1", "Allele 2", "Reads 2"):
        assert wanted in joined, f"missing column {wanted!r} in {joined}"
    assert any(h.startswith("Sequence") for h in headers)
    assert sum(h.startswith("Sequence") for h in headers) == 2, "one sequence column per allele"

    first = _profile_table(html).find(".//tbody/tr")
    assert "S-ROW" in first.text_content()


def test_sequences_are_in_their_own_scrollable_cell(tmp_path: Path) -> None:
    """A ~250 bp consensus wrapped over ten lines destroys the table; truncating
    it would hide what a sequence-resolved caller exists to show. So: scroll."""
    out = tmp_path / "r.html"
    build_report(_make_results(), RunContext(sample_name="S", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")

    boxes = _profile_table(html).find_class("seq-scroll")
    assert boxes, "sequence cells must use the scrollable container"
    assert "AATG" in boxes[0].text_content()
    css = (Path(__file__).parents[1] / "frontstr/report/static/styles.css").read_text()
    assert ".seq-scroll" in css
    block = css.split(".seq-scroll {", 1)[1].split("}", 1)[0]
    assert "overflow-x: auto" in block
    assert "white-space: nowrap" in block


def test_flags_appear_in_the_profile_table(tmp_path: Path) -> None:
    """The bug this fixes: flags existed only in the expandable per-locus cards,
    so a reviewer scanning the profile table could not see a locus was flagged."""
    out = tmp_path / "r.html"
    build_report(_flagged_results(), RunContext(sample_name="S-QC", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")
    table = _profile_table(html)

    assert "allele_imbalance" in table.text_content()
    assert table.find_class("sev-warn"), "flag chips must be coloured by severity"
    row = table.find(".//tbody/tr")
    assert "row-sev-warn" in (row.get("class") or ""), "a flagged row must be tinted"


def test_a_clean_locus_shows_no_pass_label(tmp_path: Path) -> None:
    """No aggregated PASS: a label on almost every row stops being read."""
    out = tmp_path / "r.html"
    build_report(_make_results(), RunContext(sample_name="S", panel_name="P"), out)
    table = _profile_table(out.read_text(encoding="utf-8"))
    assert "PASS" not in table.text_content()


def test_allele_balance_column_is_present(tmp_path: Path) -> None:
    out = tmp_path / "r.html"
    build_report(_flagged_results(), RunContext(sample_name="S", panel_name="P"), out)
    table = _profile_table(out.read_text(encoding="utf-8"))
    headers = " ".join(th.text_content() for th in table.findall(".//thead/tr/th"))
    assert "AB" in headers
    # 20 vs 6 reads -> 0.77, above the balanced band, so it is highlighted.
    assert table.find_class("ab-uneven")


def test_rows_are_filterable_by_flag(tmp_path: Path) -> None:
    """Both paths: the select, and typing a flag code into the search box."""
    out = tmp_path / "r.html"
    build_report(_flagged_results(), RunContext(sample_name="S", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")
    assert 'id="profile-flagged"' in html
    row = _profile_table(html).find(".//tbody/tr")
    assert row.get("data-flagged") == "1"
    assert "allele_imbalance" in (row.get("data-search") or "")
