"""Tests for the per-locus narrative trace.

The trace exists so a genotype can be followed back to the reads. These tests
guard the two properties that makes possible: that no read goes unaccounted for,
and that every stage's numbers actually appear in the rendered text.
"""

from __future__ import annotations

import pytest

from frontstr.evidence.pileup import PileupCounts, RejectReason
from frontstr.trace import (
    BinTrace,
    ClusterTrace,
    LocusTrace,
    render_locus,
    render_run_summary,
)


def _counts(fetched: int = 41, kept: int = 35, **rejected: int) -> PileupCounts:
    c = PileupCounts(fetched=fetched, kept=kept)
    for name, n in rejected.items():
        c.rejected[RejectReason[name.upper()]] = n
    return c


def _cluster(**kw: object) -> ClusterTrace:
    base: dict[str, object] = {
        "index": 0,
        "n_reads": 31,
        "fraction": 0.886,
        "length_bp": 266,
        "consensus_method": "poa_spoa",
        "n_hp1": 8,
        "n_hp2": 12,
        "n_untagged": 11,
        "number_label": "14",
        "number_method": "strnaming",
        "naming_status": "ok",
        "strnaming_name": "CE14_TATC[2]TGTC[2]TATC[10]",
        "isfg": "[TCTA]14",
        "expected_stutter": 0.0,
        "status": "allele",
        "n_reads_absorbed": 0,
        "hp_rescued": False,
        "called": True,
    }
    base.update(kw)
    return ClusterTrace(**base)  # type: ignore[arg-type]


def _trace(**kw: object) -> LocusTrace:
    base: dict[str, object] = {
        "marker": "D3S1358",
        "chrom": "chr3",
        "start": 45540634,
        "end": 45540907,
        "motif": "TCTA,TCTG",
        "period": -1,
        "strand": "+",
    }
    base.update(kw)
    return LocusTrace(**base)  # type: ignore[arg-type]


class TestReadFunnel:
    """The point of the funnel is that the arithmetic closes."""

    def test_every_fetched_read_is_accounted_for(self) -> None:
        c = _counts(fetched=41, kept=35, not_primary=3, low_mapq=1, left_flank_short=2)
        assert c.kept + c.n_rejected == c.fetched

    def test_rejection_reasons_are_rendered_with_counts(self) -> None:
        t = _trace(counts=_counts(fetched=41, kept=35, not_primary=3, left_flank_short=3))
        out = render_locus(t)
        assert "Reads fetched around the window" in out
        assert "not a primary alignment" in out
        assert "does not reach the left flank anchor" in out
        assert "35" in out

    def test_reasons_are_ordered_by_frequency(self) -> None:
        c = _counts(not_primary=1, left_flank_short=5)
        assert next(r for r, _ in c.reasons()) is RejectReason.LEFT_FLANK_SHORT

    def test_a_long_label_still_separates_from_its_value(self) -> None:
        """Regression: ljust alone ran the longest reason into its count."""
        t = _trace(counts=_counts(left_flank_short=2))
        line = next(x for x in render_locus(t).splitlines() if "left flank anchor" in x)
        assert line.rstrip().endswith("2")
        assert "anchor2" not in line

    def test_zero_coverage_stops_the_narrative_early(self) -> None:
        out = render_locus(_trace(counts=_counts(fetched=4, kept=0, low_mapq=4)))
        assert "no_data" in out
        assert "Binned by" not in out


class TestNarrativeContent:
    def test_each_stage_appears_in_order(self) -> None:
        t = _trace(
            counts=_counts(),
            bins=[BinTrace(56, 32), BinTrace(44, 1)],
            consensus_backend="poa_spoa",
            clusters=[_cluster()],
            called_labels=["14"],
            call_rule="homozygous",
        )
        out = render_locus(t)
        order = [
            out.index("Reads fetched"),
            out.index("Binned by"),
            out.index("Clustered within bins"),
            out.index("Consensus per cluster"),
            out.index("Candidates"),
            out.index("Genotype"),
        ]
        assert order == sorted(order), "the trace must follow the pipeline's order"

    def test_called_candidates_are_marked(self) -> None:
        t = _trace(
            counts=_counts(),
            clusters=[_cluster(index=0, called=True), _cluster(index=1, called=False)],
        )
        lines = [x for x in render_locus(t).splitlines() if x.strip().startswith(("*", "#"))]
        assert lines[0].strip().startswith("*")
        assert not lines[1].strip().startswith("*")

    @pytest.mark.parametrize(
        ("status", "expected"),
        [
            ("noise", "analytical threshold"),
            ("artefact", "calling threshold"),
            ("stutter", "expected as stutter"),
            ("hp_phantom", "same haplotype"),
        ],
    )
    def test_a_discarded_candidate_says_why(self, status: str, expected: str) -> None:
        """ "artefact" alone is jargon; the reader needs the rule that fired."""
        t = _trace(counts=_counts(), clusters=[_cluster(status=status, called=False)])
        assert expected in render_locus(t)

    def test_phasing_rescue_is_explained(self) -> None:
        t = _trace(counts=_counts(), clusters=[_cluster(hp_rescued=True)])
        assert "called on phasing" in render_locus(t)

    def test_unpolished_consensus_is_called_out(self) -> None:
        t = _trace(counts=_counts(), clusters=[_cluster(consensus_method="single")])
        assert "not polished by POA" in render_locus(t)

    def test_legacy_naming_path_names_the_reason(self) -> None:
        t = _trace(
            counts=_counts(),
            clusters=[
                _cluster(number_method="period_ce", naming_status="no_range", strnaming_name="")
            ],
        )
        out = render_locus(t)
        assert "period_ce" in out and "no_range" in out

    def test_a_very_long_name_is_truncated_not_wrapped(self) -> None:
        t = _trace(counts=_counts(), clusters=[_cluster(strnaming_name="CE9_" + "A" * 300)])
        assert max(len(x) for x in render_locus(t).splitlines()) < 100

    def test_markers_that_bypass_the_str_path_say_so(self) -> None:
        t = _trace(
            marker="AMEL",
            note="Sex typing, not a tandem repeat: counts reads at AMELX and AMELY.",
            called_labels=["X", "Y"],
            call_rule="heterozygous",
        )
        out = render_locus(t)
        assert "not a tandem repeat" in out
        assert "Binned by" not in out
        assert "X, Y" in out


class TestRunSummary:
    def test_totals_close_over_the_run(self) -> None:
        traces = [
            _trace(counts=_counts(fetched=40, kept=35), called_labels=["14"]),
            _trace(marker="vWA", counts=_counts(fetched=45, kept=38), called_labels=["14", "16"]),
        ]
        out = render_run_summary(traces)
        assert "85 → 73" in out
        assert "2/2" in out

    def test_loci_without_a_genotype_are_named(self) -> None:
        traces = [_trace(counts=_counts(), called_labels=[]), _trace(marker="TPOX")]
        out = render_run_summary(traces)
        assert "No genotype" in out
        assert "D3S1358" in out and "TPOX" in out

    def test_legacy_path_calls_are_surfaced(self) -> None:
        t = _trace(
            counts=_counts(),
            clusters=[_cluster(number_method="period_ce", naming_status="no_range")],
            called_labels=["14"],
        )
        out = render_run_summary([t])
        assert "legacy CE path" in out
        assert "no_range" in out

    def test_empty_run_does_not_crash(self) -> None:
        assert "No loci" in render_run_summary([])
