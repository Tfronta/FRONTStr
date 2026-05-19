"""Tests for :mod:`frontstr.report.ngs_display`."""

from __future__ import annotations

from frontstr.report.ngs_display import (
    build_ngs_panel,
    build_strhub_projection,
    highlight_repeat_spans,
)


def test_highlight_repeat_spans_wraps_longest_run() -> None:
    seq = "AAA" + "TAGA" * 4 + "CCC"
    html = highlight_repeat_spans(seq, "TAGA")
    assert "repeat-highlight" in html
    assert "seq-flank" in html
    assert '<span class="seq-flank">AAA</span>' in html
    assert '<span class="repeat-highlight">' + "TAGA" * 4 + "</span>" in html
    assert '<span class="seq-flank">CCC</span>' in html


def test_highlight_repeat_spans_all_flank_when_no_motif_match() -> None:
    html = highlight_repeat_spans("ACGTACGT", "TAGA")
    assert "seq-flank" in html
    assert "repeat-highlight" not in html
def test_build_ngs_panel_isoallele_order_and_stack() -> None:
    marker = {
        "marker_name": "vWA",
        "total_reads": 42,
        "system": {"motif": "TAGA", "period": 4},
        "alleles_called": [
            {
                "cluster_index": 0,
                "ce": 17.0,
                "consensus": "AAA" + "TAGA" * 17 + "CCCC",
                "isfg": "[TAGA]17",
                "n_reads_total": 19,
                "fraction": 19 / 42,
                "status": "allele",
                "length_bp": 75,
            },
            {
                "cluster_index": 1,
                "ce": 17.0,
                "consensus": "TTT" + "TAGA" * 17 + "CCCC",
                "isfg": "[TAGA]17",
                "n_reads_total": 23,
                "fraction": 23 / 42,
                "status": "allele",
                "length_bp": 75,
            },
        ],
        "alleles": [],
    }
    panel = build_ngs_panel(marker)
    assert len(panel["rows"]) == 2
    assert panel["rows"][0]["is_isoallele"] is True
    assert panel["rows"][1]["is_isoallele"] is False
    assert panel["rows"][0]["coverage_reads"] == 19
    assert panel["rows"][1]["coverage_reads"] == 23

    grp = panel["chart_groups"][0]
    assert grp["repeat_group"] == 17
    assert len(grp["segments"]) == 2
    assert grp["segments"][0]["coverage_reads"] == 19
    assert grp["segments"][1]["coverage_reads"] == 23


def test_build_strhub_projection_roundtrip_keys() -> None:
    marker = {
        "marker_name": "TH01",
        "total_reads": 100,
        "system": {"motif": "AATG", "period": 4},
        "alleles_called": [
            {
                "cluster_index": 0,
                "ce": 9.0,
                "consensus": "AATG" * 9,
                "isfg": "[AATG]9",
                "n_reads_total": 60,
                "fraction": 0.6,
                "status": "allele",
                "length_bp": 36,
            },
        ],
        "alleles": [],
    }
    marker["ngs_panel"] = build_ngs_panel(marker)
    hub = build_strhub_projection({"meta": {"app_version": "0", "sample_name": "S"}, "results": [marker]})
    assert hub["schema"] == "strhub.ngs_panel/v1"
    assert hub["sample"] == "S"
    assert hub["markers"][0]["locus"] == "TH01"
    row0 = hub["markers"][0]["rows"][0]
    assert row0["allele"] == "9"
    assert row0["allele_iso"] is False
    assert row0["coverage_reads"] == 60
