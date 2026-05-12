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
from frontstr.panel.models import System
from frontstr.report import RunContext, build_report


def _allele(idx: int, ce: float, cov: int, status: AlleleStatus,
            consensus_motif: str = "AATG") -> Allele:
    return Allele(
        cluster_index=idx, consensus=consensus_motif * int(ce),
        length_bp=len(consensus_motif) * int(ce),
        n_reads_total=cov, n_reads_hp1=cov // 2, n_reads_hp2=cov - cov // 2,
        n_reads_hp_none=0, n_forward=cov, n_reverse=0,
        mean_qual=30.0, ce=ce, isfg=f"[{consensus_motif}]{int(ce)}",
        bp_diff=0, is_deletion=False, status=status,
    )


def _make_results() -> list[MarkerResult]:
    th01_system = System(
        name="TH01", chromosome="chr11", ref_start=2_171_000, ref_end=2_171_050,
        motif="AATG", period=4,
    )
    tpox_system = System(
        name="TPOX", chromosome="chr2", ref_start=1_489_651, ref_end=1_489_684,
        motif="AATG", period=4,
        allow_triallelic=True, tri_balanced_thr=0.5,
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
            marker_name="TH01", system=th01_system,
            alleles=th01_alleles, alleles_called=th01_alleles[:2],
            call_rule=CallRule.HETEROZYGOUS, tri_type=TriType.NONE,
            total_reads=sum(a.n_reads_total for a in th01_alleles),
        ),
        MarkerResult(
            marker_name="TPOX", system=tpox_system,
            alleles=tpox_alleles, alleles_called=tpox_alleles,
            call_rule=CallRule.TRIALLELIC_TYPE_II, tri_type=TriType.TYPE_II_BALANCED,
            total_reads=sum(a.n_reads_total for a in tpox_alleles),
        ),
    ]


def test_build_report_writes_valid_html(tmp_path: Path) -> None:
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(sample_name="S001", panel_name="forensic-panel", panel_version="0.1",
                   operator="J. Diaz"),
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

    # Required sections present
    for sid in ("cover", "profile", "qc", "loci", "raw"):
        assert doc.get_element_by_id(sid) is not None

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
    sha_in_footer = re.search(r"id=\"report-hash-footer\">([0-9a-f]{64})<", html)
    assert sha_in_footer is not None


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
    assert html.count("class=\"chip ok\"") > 0
    assert "triallelic_type_II" in html
