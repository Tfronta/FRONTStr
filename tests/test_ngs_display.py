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
    hub = build_strhub_projection(
        {"meta": {"app_version": "0", "sample_name": "S"}, "results": [marker]}
    )
    assert hub["schema"] == "strhub.ngs_panel/v1"
    assert hub["sample"] == "S"
    assert hub["markers"][0]["locus"] == "TH01"
    row0 = hub["markers"][0]["rows"][0]
    assert row0["allele"] == "9"
    assert row0["allele_iso"] is False
    assert row0["coverage_reads"] == 60


def test_chart_group_is_labelled_with_the_allele_not_the_rounded_group() -> None:
    """A 9.3 microvariant must not be ticked "9" on the chart.

    ``repeat_group`` is round(CE). It exists so isoalleles stack on one bar and
    is not an allele designation. Using it as the tick label put "9" under the
    TH01 microvariant, contradicting the table directly above the chart and
    erasing the very thing that makes the locus interesting.
    """
    marker = {
        "marker_name": "TH01",
        "total_reads": 17,
        "system": {"motif": "AATG", "period": 4},
        "alleles_called": [
            {
                "cluster_index": 0,
                "ce": 9.3,
                "consensus": "AAA" + "AATG" * 9 + "ATG" + "CCC",
                "isfg": "[AATG]9 ATG",
                "n_reads_total": 10,
                "fraction": 10 / 17,
                "status": "allele",
                "length_bp": 39,
            },
            {
                "cluster_index": 1,
                "ce": 7.0,
                "consensus": "AAA" + "AATG" * 7 + "CCC",
                "isfg": "[AATG]7",
                "n_reads_total": 7,
                "fraction": 7 / 17,
                "status": "allele",
                "length_bp": 28,
            },
        ],
    }
    panel = build_ngs_panel(marker)
    labels = {g["label"] for g in panel["chart_groups"]}
    assert labels == {"9.3", "7"}
    assert "9" not in labels  # the rounded group, which is what it used to show

    # And the tick agrees with the table row for the same allele.
    displays = {r["allele_display"] for r in panel["rows"]}
    assert labels <= displays


def test_chart_group_label_joins_distinct_designations() -> None:
    """One bar may hold more than one designation; the tick names them all."""
    marker = {
        "marker_name": "D3S1358",
        "total_reads": 30,
        "system": {"motif": "TCTA", "period": 4},
        "alleles_called": [
            {
                "cluster_index": 0,
                "ce": 15.0,
                "consensus": "AAA" + "TCTA" * 15 + "CCC",
                "isfg": "[TCTA]15",
                "n_reads_total": 18,
                "fraction": 0.6,
                "status": "allele",
                "length_bp": 60,
            },
            {
                "cluster_index": 1,
                "ce": 15.0,
                "consensus": "AAA" + "TCTA" * 14 + "TCTG" + "CCC",
                "isfg": "[TCTA]14 TCTG",
                "n_reads_total": 12,
                "fraction": 0.4,
                "status": "allele",
                "length_bp": 60,
            },
        ],
    }
    groups = build_ngs_panel(marker)["chart_groups"]
    assert len(groups) == 1, "same repeat group must stay one bar"
    assert len(groups[0]["segments"]) == 2, "isoalleles stack on that bar"
    assert groups[0]["label"] == "15", "one designation, stated once, not '15/15'"
