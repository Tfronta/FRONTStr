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
    assert "Read clusters" in html
    assert "sequence-scroll" in html


def test_report_does_not_claim_a_forensic_audit_was_performed(tmp_path: Path) -> None:
    """The per-locus cluster table was headed "Forensic audit".

    It is a table of read clusters. A forensic audit is a review carried out by
    people against a standard, and the README states in as many words that the
    caller is not laboratory-validated and that no forensic partner has signed
    off. A heading asserting one had happened is the kind of claim that gets a
    tool disqualified when someone reads it literally. "Raw / audit" stays: an
    audit *record* is the trail that would let one be performed, which is what
    the hashes and the parameter list are.
    """
    out = tmp_path / "report.html"
    build_report(_make_results(), RunContext(sample_name="S", panel_name="P"), out)
    assert "forensic audit" not in out.read_text(encoding="utf-8").lower()


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
    # Forensic nomenclature: AD is allelic depth, DP is the depth behind the
    # call. "Reads 1" said the same thing in words nobody reports genotypes in.
    for wanted in ("Sample", "Marker", "Allele 1", "Allele 2", "AD"):
        assert wanted in joined, f"missing column {wanted!r} in {joined}"
    assert "Reads" not in joined, "the old label must be gone, not merely joined by a new one"
    # No locus total. Depth is reported per allele and only per allele: that is
    # the evidence behind each allele separately, and it is what clustering by
    # sequence buys. A summed column also invites the discarded reads to be read
    # as if they backed the genotype.
    for gone in ("Cov", "DP", "ΣAD", "Total"):
        assert not any(h.strip().startswith(gone) for h in headers), (
            f"{gone} is a locus total; depth belongs per allele"
        )
    # The bracketed ISFG string, not the raw consensus. This table is the
    # profile a reader compares against another profile, and a wall of bases
    # cannot be compared by eye; the full consensus is in Sequences below,
    # where the row is per allele and there is room for it.
    assert sum(h.startswith("ISFG") for h in headers) == 2, "one ISFG column per allele"
    assert not any(h.startswith("Sequence") for h in headers), (
        "the raw consensus belongs in the Sequences section, not the profile"
    )

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

    # The chip is the abbreviation, so the QC column stays narrow enough that
    # Allele 2 and its read count do not fall off the right edge. The full code
    # rides in the tooltip and in the legend under the table.
    assert "AI" in table.text_content(), "the abbreviated code must be on the chip"
    assert "allele_imbalance" not in table.text_content(), (
        "the full code belongs in the tooltip, not in the cell"
    )
    assert 'title="allele_imbalance' in html, "the tooltip must carry the full code"
    assert table.find_class("sev-warn"), "flag chips must be coloured by severity"
    row = table.find(".//tbody/tr")
    assert "row-sev-warn" in (row.get("class") or ""), "a flagged row must be tinted"


def test_qc_legend_lists_only_the_codes_this_run_raised(tmp_path: Path) -> None:
    """A legend of all fifteen codes is a wall nobody reads.

    It also has to expand the abbreviations actually shown, or the narrow chips
    become unreadable jargon.
    """
    out = tmp_path / "r.html"
    build_report(_flagged_results(), RunContext(sample_name="S-QC", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")

    legend = lxml.html.fromstring(html).find_class("qc-legend")
    assert legend, "no QC legend was rendered"
    text = legend[0].text_content()
    assert "allele_imbalance" in text, "the raised code must be expanded"
    assert "phase_block_split" not in text, "a code that never fired must not be listed"


def test_every_flag_code_has_an_abbreviation() -> None:
    """A new code without one would render a blank chip, silently."""
    from frontstr.interp.models import FlagCode

    shorts = {code: code.short for code in FlagCode}

    assert all(shorts.values()), "every code needs a non-empty abbreviation"
    assert len(set(shorts.values())) == len(shorts), f"abbreviations collide: {shorts}"
    assert all(len(s) <= 3 for s in shorts.values()), "abbreviations must stay narrow"


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


# ---------------------------------------------------------------------------
# One repeat string everywhere, and run provenance
# ---------------------------------------------------------------------------


def test_report_shows_the_strnaming_string_not_the_window_scan(tmp_path: Path) -> None:
    """The bug: --trace printed CE9.3_TGAA[6]… while the HTML printed the
    legacy full-window scan, a hundred lowercase flank bases before the
    brackets. Two strings for one allele, in a forensic report."""
    alleles = [_allele(0, 9.0, 20, AlleleStatus.ALLELE)]
    alleles[0].strnaming_name = "CE9_TGAA[9]"
    alleles[0].strnaming_ce = 9.0
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
        call_rule=CallRule.HOMOZYGOUS,
        tri_type=TriType.NONE,
        total_reads=20,
    )
    out = tmp_path / "r.html"
    build_report([r], RunContext(sample_name="S", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")
    assert "CE9_TGAA[9]" in html, "the STRNaming name must reach the report"


def test_markers_without_a_strnaming_range_keep_the_legacy_string(tmp_path: Path) -> None:
    """DYS393 and AMEL have no reporting range; they must still show something."""
    out = tmp_path / "r.html"
    build_report(_make_results(), RunContext(sample_name="S", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")
    assert "[AATG]9" in html, "the bracket-scan fallback must still render"


def test_report_embeds_the_panel_bed(tmp_path: Path) -> None:
    """Naming the loci is not the same as showing the intervals they came from."""
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(
            sample_name="S",
            panel_name="P",
            panel_bed=["chr11\t2171000\t2171050\tAATG\tTH01"],
        ),
        out,
    )
    html = out.read_text(encoding="utf-8")
    assert "Panel windows (BED)" in html
    assert "2171000" in html and "2171050" in html


def test_report_shows_the_command_and_every_effective_parameter(tmp_path: Path) -> None:
    """argv alone is misleading: the values that decide a call are the defaults."""
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(
            sample_name="S",
            panel_name="P",
            pipeline_argv=["frontstr", "export", "--bam", "x.bam"],
            effective_params={"min_mapq": 20, "calling_thresh": 0.1},
        ),
        out,
    )
    html = out.read_text(encoding="utf-8")
    assert "frontstr export --bam x.bam" in html
    assert "Parameters in force" in html
    assert "min_mapq" in html and "calling_thresh" in html


def test_provenance_sections_are_absent_when_unknown(tmp_path: Path) -> None:
    """A library caller supplies none of this; the report must not fake it."""
    out = tmp_path / "r.html"
    build_report(_make_results(), RunContext(sample_name="S", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")
    assert "Panel windows (BED)" not in html
    assert "Parameters in force" not in html


def test_profile_reports_depth_for_every_called_allele(tmp_path: Path) -> None:
    """Every called allele carries its own read count, and nothing sums them.

    Per-allele depth is the differentiator: a caller that only sizes the repeat
    cannot say how many reads stand behind each allele. Regression guard for the
    column that used to total them, which read as if reads the caller had
    discarded were part of the evidence.
    """
    out = tmp_path / "r.html"
    build_report(_make_results(), RunContext(sample_name="S-COV", panel_name="P"), out)
    table = _profile_table(out.read_text(encoding="utf-8"))

    headers = [th.text_content().strip() for th in table.findall(".//thead/tr/th")]
    ad_cols = [i for i, h in enumerate(headers) if h.startswith("AD")]
    assert len(ad_cols) >= 2, "one AD column per allele slot"

    checked = 0
    for tr in table.findall(".//tbody/tr"):
        tds = tr.findall("td")
        for i in ad_cols:
            # The cell holds the depth plus a .strand-split child; read the
            # depth off the cell's own text so the two cannot be conflated.
            depth = (tds[i].text or "").strip()
            if depth.isdigit():
                assert int(depth) > 0, "a called allele with zero depth is not a call"
                checked += 1
    assert checked, "no per-allele depth was rendered at all"


def test_strand_balance_is_its_own_table_not_a_column_on_the_profile(tmp_path: Path) -> None:
    """Strand evidence gets a view of its own, with a header that names it.

    It was briefly a sixth quantity inside the profile table, which already
    carries allele number, depth, balance, QC and ISFG per row. Too much at
    once, and unlabelled: the reader could not tell the pair apart from the
    depth it sat under. ``STRAND_BIAS`` cannot speak below
    ``strand_bias_min_reads``, and the thin alleles it cannot reach are the
    ones most likely to be artefacts, so the numbers have to be somewhere.
    """
    out = tmp_path / "r.html"
    build_report(_make_results(), RunContext(sample_name="S-STR", panel_name="P"), out)
    html = out.read_text(encoding="utf-8")

    # Not on the profile table.
    assert not _profile_table(html).findall('.//span[@class="strand-split"]')

    root = lxml.html.fromstring(html)
    heads = [h.text_content().strip() for h in root.findall(".//h3")]
    assert "Strand balance" in heads, "the strand table has no heading of its own"

    tables = root.findall('.//table[@data-sortable-table]')
    assert tables, "the strand table is not sortable"
    headers = [th.text_content().strip() for th in tables[0].findall(".//thead/tr/th")]
    assert headers[:5] == ["Marker", "Allele", "Reads", "Fwd", "Rev"]
    # The locus is carried beside the allele: it is the baseline that makes the
    # allele's split readable, and a locus can be skewed on its own account.
    assert sum("whole locus" in h for h in headers) == 1
    assert tables[0].findall(".//tbody/tr"), "no called alleles listed"


def test_no_way_back_to_the_cohort_without_one(tmp_path: Path) -> None:
    """A standalone run has nowhere to go back to, so it gets no link."""
    out = tmp_path / "r.html"
    build_report(_make_results(), RunContext(sample_name="S", panel_name="P"), out)

    # Against the DOM, not the raw text: the class name is also in the inlined
    # stylesheet, which is always present whether or not the link is rendered.
    doc = lxml.html.fromstring(out.read_text(encoding="utf-8"))
    assert not doc.find_class("nav-back"), "a link to a cohort view that does not exist"


def test_a_sample_from_a_cohort_can_get_back_to_it(tmp_path: Path) -> None:
    """Reached from the cohort table, the report has to lead back.

    Without it the only way out was the browser's back button, which is not a
    control anyone should have to discover.
    """
    out = tmp_path / "r.html"
    build_report(
        _make_results(),
        RunContext(sample_name="S", panel_name="P", cohort_href="../cohort.html"),
        out,
    )
    doc = lxml.html.fromstring(out.read_text(encoding="utf-8"))

    back = doc.find_class("nav-back")
    assert back, "no link back to the cohort"
    assert back[0].get("href") == "../cohort.html"
