"""End-to-end regression against a real ONT sample with a known profile.

Every other test in this suite runs on synthetic fixtures. This one runs the
whole pipeline — pileup, clustering, POA consensus, haplotype suppression,
stutter, ISFG, nomenclature — over a real ONT R10 BAM and checks the genotypes
against a curated reference profile. Without it, a change can keep 300 unit
tests green while silently changing what FRONTStr calls.

The BAM slices are not versioned (they are ~200 MB), so these tests skip when
the data is absent. To run them::

    pytest tests/test_regression_hg00113.py -m integration

Sample: HG00113, 1000 Genomes GBR, male. Slice of
``HG00113-ONT-hg38-R10-LSK114-dorado090_sup_5mCG_5hmCG_v500.phased.bam``,
±10 kb around each panel locus.

Reference genotypes come from the ``Illumina`` sheet of
``1000GEN-ONT-Merged-Compar.xlsx`` (HipSTR on the matched Illumina data). Its
``Concordancia-3-tecnologias`` sheet adds LongTR and STRspy side by side and is
what adjudicates a disagreement. Two caveats that matter when reading a failure
here:

- **D21S11 is not called by HipSTR at all** (``NA`` for every sample), so its
  29/31 rests on LongTR and STRspy only — both long-read, i.e. caller-vs-caller
  rather than an orthogonal method.
- **The sex markers have no external reference whatsoever.** DYS391 is ``NA``
  and DYS393/DXS7132/DXS8378 are absent from the workbook, so those four
  expectations are FRONTStr's own output. Do not treat them as validated.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.interp.models import CallRule, FlagCode, MarkerResult
from frontstr.interp.profile import interpret_run
from frontstr.panel.loader import load_panel

SLICE_DIR = Path(__file__).parent / "data" / "ont_slices"
PANEL_PATH = Path(__file__).parents[1] / "examples" / "panels" / "codis_20_grch38.yaml"

pytestmark = [pytest.mark.integration, pytest.mark.slow]


#: Expected allele-number labels per marker, as ``Allele.number_label`` renders
#: them. Sets, because allele order within a genotype is not meaningful.
#:
#: This table used to carry a documented deviation at vWA (13/17) and D21S11
#: (35/33), where FRONTStr's bracket repeat count did not map onto the kit
#: designation and no single ``corr_value`` could reconcile it. Adopting
#: STRNaming as the source of the allele number (:mod:`frontstr.interp.naming`)
#: resolved both, so **every marker here is now the kit/Book1 value**. DXS7132
#: moved 15 → 14 in the same change: its old ``corr_value`` had been calibrated
#: against FRONTStr's own output rather than an external truth, and STRNaming
#: agrees with the GRCh38 reference designation.
EXPECTED_PROFILE: dict[str, set[str]] = {
    "D3S1358": {"14"},
    "vWA": {"14", "16"},
    "FGA": {"24", "21"},
    "D8S1179": {"13", "10"},
    "D21S11": {"29", "31"},
    "D18S51": {"14", "15"},
    "D5S818": {"11", "13"},
    "D13S317": {"11", "9"},
    "D7S820": {"8", "10"},
    "D16S539": {"11", "8"},
    "TH01": {"9.3", "7"},
    "TPOX": {"9", "11"},
    "CSF1PO": {"10", "12"},
    # 17 is carried by 5 reads against the major allele's 17 — a PHR of 0.29,
    # under the 0.4 floor. Called heterozygous because the two clusters are
    # 100% HP1 and 100% HP2 respectively; see the HP-rescue test below.
    "D2S1338": {"20", "17"},
    "D19S433": {"14", "12"},
    "D10S1248": {"13", "15"},
    "D1S1656": {"13", "11"},
    "D2S441": {"13", "10"},
    "D12S391": {"21", "20"},
    "D22S1045": {"16"},
    "AMEL": {"X", "Y"},
    # Hemizygous in a male; reported as homozygous today (see ROADMAP).
    "DYS391": {"10"},
    "DYS393": {"14"},
    "DXS7132": {"14"},
    "DXS8378": {"10"},
}

#: Markers STRNaming defines no reported range for, which therefore keep the
#: legacy CE path. AMEL is not an STR; DYS393 is simply absent from STRNaming's
#: ``ranges_uas-frr.txt``.
NO_STRNAMING_RANGE = {"AMEL", "DYS393"}

#: Single-source 1000 Genomes individuals. None of these can be a true mixture,
#: so any ``mixture_suspected`` here is a false positive by construction.
SINGLE_SOURCE_SAMPLES = (
    "HG00097",
    "HG00113",
    "HG00154",
    "HG00263",
    "GM19038",
)


def _slice_path(sample: str) -> Path:
    return SLICE_DIR / f"{sample}.codis.bam"


def _require(sample: str) -> Path:
    path = _slice_path(sample)
    if not path.exists():
        pytest.skip(f"ONT slice not available: {path} (not versioned — see docstring)")
    return path


@pytest.fixture(scope="module")
def hg00113_results() -> list[MarkerResult]:
    bam = _require("HG00113")
    return interpret_run(bam=bam, panel=load_panel(PANEL_PATH))


def _called_labels(result: MarkerResult) -> set[str]:
    return {a.number_label for a in result.alleles_called}


def test_every_panel_marker_produces_a_result(hg00113_results: list[MarkerResult]) -> None:
    assert {r.marker_name for r in hg00113_results} == set(EXPECTED_PROFILE)


@pytest.mark.parametrize("marker", sorted(EXPECTED_PROFILE))
def test_marker_genotype_matches_reference(
    hg00113_results: list[MarkerResult], marker: str
) -> None:
    result = next(r for r in hg00113_results if r.marker_name == marker)
    assert _called_labels(result) == EXPECTED_PROFILE[marker], (
        f"{marker}: called {sorted(_called_labels(result))}, "
        f"expected {sorted(EXPECTED_PROFILE[marker])}"
    )


def test_no_locus_drops_out(hg00113_results: list[MarkerResult]) -> None:
    """Every marker in the panel is covered in this slice; none may go no_data."""
    dropouts = [r.marker_name for r in hg00113_results if r.call_rule == CallRule.NO_DATA]
    assert dropouts == []


def test_no_false_mixture_flags(hg00113_results: list[MarkerResult]) -> None:
    """HG00113 is a single-source individual — see known-bug #6."""
    flagged = [
        r.marker_name
        for r in hg00113_results
        if any(f.code == FlagCode.MIXTURE_SUSPECTED for f in r.flags)
    ]
    assert flagged == []


def test_consensus_is_polished(hg00113_results: list[MarkerResult]) -> None:
    """A POA backend must be active: an unpolished consensus fakes microvariants.

    Measured on this slice set, the mode fallback produced 4 false microvariants
    in 202 called alleles. If this fails, install ``frontstr[poa]``.
    """
    unpolished = [
        r.marker_name
        for r in hg00113_results
        if any(f.code == FlagCode.CONSENSUS_FALLBACK for f in r.flags)
    ]
    assert unpolished == []


def test_th01_microvariant_survives_clustering(hg00113_results: list[MarkerResult]) -> None:
    """TH01 9.3 is the headline NGS-over-CE claim; it must not be merged into 9 or 10.

    Repeat-core binning is what makes this safe: 9 and 9.3 have the same repeat
    unit count (9) and differ only in core length (36 vs 39 bp).
    """
    th01 = next(r for r in hg00113_results if r.marker_name == "TH01")
    nine_three = next(a for a in th01.alleles_called if a.number_label == "9.3")
    assert nine_three.isfg.count("AATG") >= 2, "9.3 must keep its interrupted structure"
    assert "ATG " in nine_three.isfg or nine_three.isfg.endswith("ATG")


def test_allele_numbers_come_from_strnaming(hg00113_results: list[MarkerResult]) -> None:
    """Every marker with a defined reporting range must be named by STRNaming.

    Guards the wiring end to end: a missing slice cache, a panel window too
    narrow to contain a range, or an anchor failure would silently drop that
    marker back to the legacy bracket count — which is exactly the arithmetic
    that was wrong at six loci.
    """
    fell_back = {
        r.marker_name: sorted({a.strnaming_status for a in r.alleles_called})
        for r in hg00113_results
        if r.marker_name not in NO_STRNAMING_RANGE
        and any(a.number_method != "strnaming" for a in r.alleles_called)
    }
    assert fell_back == {}, f"markers silently on the legacy CE path: {fell_back}"


def test_markers_without_a_range_keep_the_legacy_path(
    hg00113_results: list[MarkerResult],
) -> None:
    """The fallback must stay working — it is the only path for these markers."""
    dys393 = next(r for r in hg00113_results if r.marker_name == "DYS393")
    assert [a.strnaming_status for a in dys393.alleles_called] == ["no_range"]
    assert {a.number_method for a in dys393.alleles_called} == {"period_ce"}
    assert {a.number_label for a in dys393.alleles_called} == EXPECTED_PROFILE["DYS393"]


def test_imbalanced_heterozygote_is_rescued_by_phasing(
    hg00113_results: list[MarkerResult],
) -> None:
    """D2S1338 is 17/20 in Illumina, LongTR and STRspy; FRONTStr used to say 20.

    The minor allele has 5 reads against 17 — a peak-height ratio of 0.29,
    under the 0.4 het floor inherited from capillary electrophoresis. The two
    clusters are 100% HP1 and 100% HP2, which is direct evidence of two alleles,
    so the ratio must not overrule it. Guards against the rule silently
    regressing to a false homozygote — a false exclusion, the costliest error
    this caller can make.
    """
    result = next(r for r in hg00113_results if r.marker_name == "D2S1338")
    assert result.call_rule == CallRule.HETEROZYGOUS
    rescued = [a for a in result.alleles_called if a.hp_rescued]
    assert len(rescued) == 1, "the 17 allele should be called on phasing"
    assert rescued[0].number_label == "17"
    assert any(f.code == FlagCode.HP_RESCUED_HET for f in result.flags), (
        "a call that rests on phasing rather than peak balance must be flagged"
    )


def test_phasing_rescue_does_not_invent_heterozygotes(
    hg00113_results: list[MarkerResult],
) -> None:
    """The rescue must stay rare and targeted, not loosen calling generally.

    It fires at exactly one locus in this sample. If a change makes it fire
    broadly, it has stopped being an appeal to phasing evidence and become a
    lowered threshold.
    """
    rescued = [
        r.marker_name for r in hg00113_results if any(a.hp_rescued for a in r.alleles_called)
    ]
    assert rescued == ["D2S1338"]


def test_resolved_markers_no_longer_warn_about_kit_nomenclature(
    hg00113_results: list[MarkerResult],
) -> None:
    """vWA and D21S11 now report the kit designation, so the warning must be gone.

    The panel still carries ``kit_nomenclature_note`` for them as documentation
    of the legacy fallback; the flag is what must not fire, because telling a
    reviewer not to compare a number that *is* the kit value would cause the
    false exclusion the flag exists to prevent.
    """
    for marker in ("vWA", "D21S11"):
        result = next(r for r in hg00113_results if r.marker_name == marker)
        offenders = [f for f in result.flags if f.code == FlagCode.CE_NOMENCLATURE_OFFSET]
        assert offenders == [], f"{marker} still warns: {[f.message for f in offenders]}"


@pytest.mark.parametrize("sample", SINGLE_SOURCE_SAMPLES)
def test_single_source_samples_raise_no_mixture(sample: str) -> None:
    """The regression lock for known-bug #6 across the whole slice set.

    Repeat-core binning removes the split-allele phantoms at source and
    haplotype-aware suppression catches the remainder. Before either, five loci
    across these samples raised ``mixture_suspected``.
    """
    bam = _require(sample)
    results = interpret_run(bam=bam, panel=load_panel(PANEL_PATH))
    flagged = [
        r.marker_name for r in results if any(f.code == FlagCode.MIXTURE_SUSPECTED for f in r.flags)
    ]
    assert flagged == [], f"{sample}: false mixture at {flagged}"


def test_female_samples_have_no_y_signal() -> None:
    """Y markers must return no_data in females under the default MAPQ filter.

    X/Y paralogues produce secondary alignments that only MAPQ filtering
    excludes, so this guards the filter as much as the sex typing.
    """
    bam = _require("HG00097")  # female
    results = interpret_run(bam=bam, panel=load_panel(PANEL_PATH))
    for marker in ("DYS391", "DYS393"):
        result = next(r for r in results if r.marker_name == marker)
        assert result.call_rule == CallRule.NO_DATA, f"{marker} called in a female"


def test_trace_accounts_for_every_read_at_every_locus() -> None:
    """The trace's arithmetic must close on real data, at every marker.

    A coverage number a reviewer cannot reconcile against the BAM is not
    auditable. This is the property the whole trace exists to provide, so it is
    checked end to end rather than only on synthetic fixtures.
    """
    from frontstr.trace import LocusTrace

    bam = _require("HG00113")
    traces: list[LocusTrace] = []
    interpret_run(bam=bam, panel=load_panel(PANEL_PATH), on_trace=traces.append)

    assert {t.marker for t in traces} == set(EXPECTED_PROFILE), "every marker must be traced"
    for t in traces:
        if t.counts is None:
            assert t.note, f"{t.marker}: untraced without saying why"
            continue
        assert t.counts.kept + t.counts.n_rejected == t.counts.fetched, (
            f"{t.marker}: {t.counts.fetched} fetched but "
            f"{t.counts.kept} kept + {t.counts.n_rejected} rejected"
        )


def test_trace_coverage_matches_the_reported_genotype() -> None:
    """The funnel's survivor count must be the coverage the profile reports."""
    from frontstr.trace import LocusTrace

    bam = _require("HG00113")
    traces: list[LocusTrace] = []
    results = interpret_run(bam=bam, panel=load_panel(PANEL_PATH), on_trace=traces.append)
    by_marker = {r.marker_name: r for r in results}

    for t in traces:
        if t.counts is None:
            continue
        assert t.counts.kept == by_marker[t.marker].total_reads, (
            f"{t.marker}: trace says {t.counts.kept} reads, profile says "
            f"{by_marker[t.marker].total_reads}"
        )


def test_phase_block_split_is_detected_and_declines_to_guess() -> None:
    """HP labels are local to a phase block, and 3 of 125 loci really do split.

    HG00097 D13S317 is the clear case: a 14-read cluster whose HP2 labels are
    100% pure but drawn from two different blocks. Nothing about the call
    changes — the haplotype rules simply decline — but the flag has to fire, or
    a reviewer reads those HP counts as haplotype evidence they are not.
    """
    from frontstr.interp.haplotype import dominant_hp

    bam = _require("HG00097")
    results = interpret_run(bam=bam, panel=load_panel(PANEL_PATH))
    d13 = next(r for r in results if r.marker_name == "D13S317")

    split = [a for a in d13.alleles if a.n_phase_sets > 1]
    assert split, "the two-block cluster must be recorded"
    assert all(dominant_hp(a) is None for a in split), (
        "a cluster spanning blocks must get no haplotype, however pure its HP labels"
    )
    assert any(f.code == FlagCode.PHASE_BLOCK_SPLIT for f in d13.flags)
    # The genotype is unaffected: this is a latent-evidence fix, not a call fix.
    assert _called_labels(d13) == {"11", "14"}


def test_phase_blocks_are_read_from_the_bam(hg00113_results: list[MarkerResult]) -> None:
    """If PS stopped being parsed, every guard above would silently pass."""
    with_blocks = [a for r in hg00113_results for a in r.alleles if a.phase_set is not None]
    assert with_blocks, "the 1000G ONT slices are phased and do carry PS"


def test_lowering_a_measured_threshold_marks_every_marker() -> None:
    """A laboratory may test a threshold; it may not do so invisibly.

    ``min_reads_third`` = 5 came from the known-bug #6 investigation into ONT
    basecaller phantoms. Lowering it re-admits those phantoms, so the resulting
    profile is not comparable with a default run — and six months later nobody
    remembers which run was the experiment. The flag is per marker because that
    is where the audit census, the XLSX QC sheet and the HTML row tint all look.
    """
    from frontstr.params import RunParameters

    bam = _require("HG00113")
    panel = load_panel(PANEL_PATH)

    clean = interpret_run(bam=bam, panel=panel, params=RunParameters.defaults())
    assert not any(f.code == FlagCode.NON_DEFAULT_THRESHOLD for r in clean for f in r.flags), (
        "a default run must not be marked"
    )

    tuned = interpret_run(bam=bam, panel=panel, params=RunParameters.of(min_reads_third=2))
    flagged = [r for r in tuned if any(f.code == FlagCode.NON_DEFAULT_THRESHOLD for f in r.flags)]
    assert len(flagged) == len(tuned), "every marker must carry the mark"
    assert "min_reads_third=2" in flagged[0].flags[-1].message
    assert "default 5" in flagged[0].flags[-1].message


def test_tuning_a_chosen_threshold_does_not_mark_the_run() -> None:
    """Only measured defaults mark. Ordinary tuning must stay quiet, or the
    flag fires so often it stops being read."""
    from frontstr.params import RunParameters

    results = interpret_run(
        bam=_require("HG00113"),
        panel=load_panel(PANEL_PATH),
        params=RunParameters.of(calling_thresh=0.08),
    )
    assert not any(f.code == FlagCode.NON_DEFAULT_THRESHOLD for r in results for f in r.flags)


def test_parameters_actually_reach_the_pipeline() -> None:
    """The knobs must change behaviour, not just appear in the echo."""
    from frontstr.params import RunParameters

    bam = _require("HG00113")
    panel = load_panel(PANEL_PATH)
    strict = interpret_run(bam=bam, panel=panel, params=RunParameters.of(min_mapq=60))
    default = interpret_run(bam=bam, panel=panel, params=RunParameters.defaults())
    assert sum(r.total_reads for r in strict) < sum(r.total_reads for r in default), (
        "raising min_mapq must drop reads"
    )
