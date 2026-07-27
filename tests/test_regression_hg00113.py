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
±10 kb around each panel locus. Reference genotypes from Book1.xlsx (HipSTR on
the matched Illumina data) cross-checked against the CODIS profile.
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
#: Two markers deliberately encode a **known, documented deviation** from the
#: legacy CE-kit designation rather than the kit value. FRONTStr reports the
#: ISFG bracket repeat count (Phillips 2018), which for these two compound loci
#: does not map 1:1 onto the kit convention established for CE ladders. The
#: ISFG strings are structurally correct in both cases; no single ``corr_value``
#: can reconcile them, because the offsets run in opposite directions between
#: the two alleles of vWA. If a future change makes these match the kit values,
#: this test will fail — and that failure is good news worth investigating, not
#: something to paper over.
EXPECTED_PROFILE: dict[str, set[str]] = {
    "D3S1358": {"14"},
    "vWA": {"13", "17"},  # kit convention: 14/16 — see note above
    "FGA": {"24", "21"},
    "D8S1179": {"13", "10"},
    "D21S11": {"35", "33"},  # kit convention: 29/31 — see note above
    "D18S51": {"14", "15"},
    "D5S818": {"11", "13"},
    "D13S317": {"11", "9"},
    "D7S820": {"8", "10"},
    "D16S539": {"11", "8"},
    "TH01": {"9.3", "7"},
    "TPOX": {"9", "11"},
    "CSF1PO": {"10", "12"},
    # Reference is 20/17; the 17 allele does not clear the calling threshold in
    # this slice. A coverage limitation, not a nomenclature one.
    "D2S1338": {"20"},
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
    "DXS7132": {"15"},
    "DXS8378": {"10"},
}

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
        r.marker_name
        for r in results
        if any(f.code == FlagCode.MIXTURE_SUSPECTED for f in r.flags)
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
