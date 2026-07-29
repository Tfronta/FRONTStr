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
    RunHeader,
    render_header,
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
        assert "grouped by length" not in out


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
            out.index("Step 1 — grouped by length"),
            out.index("Step 2 — split by sequence"),
            out.index("Step 3 — consensus per cluster"),
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
        assert "NOT polished by POA" in render_locus(t)

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
        name_line = next(x for x in render_locus(t).splitlines() if x.strip().startswith("name"))
        assert len(name_line) < 100

    def test_markers_that_bypass_the_str_path_say_so(self) -> None:
        t = _trace(
            marker="AMEL",
            note="Sex typing, not a tandem repeat: counts reads at AMELX and AMELY.",
            called_labels=["X", "Y"],
            call_rule="heterozygous",
        )
        out = render_locus(t)
        assert "not a tandem repeat" in out
        assert "grouped by length" not in out
        assert "X, Y" in out


class TestSequences:
    """The bases themselves — the point of a sequence-resolved caller."""

    def test_flank_core_flank_is_rendered(self) -> None:
        t = _trace(
            counts=_counts(),
            clusters=[_cluster(flank_left="…GACTCCATGGTG", core="AATG" * 7, flank_right="AGGG…")],
        )
        out = render_locus(t)
        assert "Sequences" in out
        assert "AATGAATGAATGAATGAATGAATGAATG" in out
        assert "…GACTCCATGGTG" in out

    def test_cores_are_column_aligned_across_candidates(self) -> None:
        """Reading down the column is how a 4 bp step becomes visible."""
        t = _trace(
            counts=_counts(),
            clusters=[
                _cluster(index=0, flank_left="…AAAA", core="AATG" * 9, flank_right="GGGG…"),
                _cluster(index=1, flank_left="…AAAA", core="AATG" * 7, flank_right="GGGG…"),
            ],
        )
        rows = [x for x in render_locus(t).splitlines() if "AATG" in x]
        assert len(rows) == 2
        assert rows[0].index("AATG") == rows[1].index("AATG"), "cores must start in one column"

    def test_a_called_allele_is_never_truncated(self) -> None:
        """D21S11's core runs past 180 bp; cutting it would defeat the purpose."""
        long_core = "TCTA" * 50  # 200 bp
        t = _trace(counts=_counts(), clusters=[_cluster(core=long_core, called=True)])
        assert long_core in render_locus(t)

    def test_an_uncalled_noise_candidate_is_truncated(self) -> None:
        t = _trace(counts=_counts(), clusters=[_cluster(core="A" * 300, called=False)])
        out = render_locus(t)
        assert "A" * 300 not in out
        assert "…" in out
        assert "never are" in out, "the cut must be disclosed, not silent"

    def test_a_missing_core_is_disclosed(self) -> None:
        t = _trace(
            counts=_counts(),
            clusters=[_cluster(core="ACGT" * 5, core_found=False, flank_left="", flank_right="")],
        )
        assert "no motif run found" in render_locus(t)

    def test_no_sequence_section_when_there_is_nothing_to_show(self) -> None:
        t = _trace(counts=_counts(), clusters=[_cluster(core="", flank_left="", flank_right="")])
        assert "Sequences" not in render_locus(t)

    def test_rows_have_no_trailing_whitespace(self) -> None:
        t = _trace(
            counts=_counts(),
            clusters=[
                _cluster(index=0, core="AATG" * 9, flank_right=""),
                _cluster(index=1, core="AATG" * 2, flank_right=""),
            ],
        )
        assert all(x == x.rstrip() for x in render_locus(t).splitlines())


class TestRunHeader:
    """A benchmark log has to state its own inputs, or it is useless later."""

    def test_counts_bams_and_crams_separately(self) -> None:
        out = render_header(RunHeader(inputs=["a.bam", "b.bam", "c.cram"]))
        assert "2 BAM" in out and "1 CRAM" in out

    def test_lists_each_input_path(self) -> None:
        out = render_header(RunHeader(inputs=["one.bam", "two.bam"]))
        assert "one.bam" in out and "two.bam" in out

    def test_states_every_threshold_in_force(self) -> None:
        out = render_header(
            RunHeader(
                inputs=["x.bam"],
                min_mapq=25,
                analytical_thresh=0.03,
                calling_thresh=0.12,
                identity_threshold=0.95,
            )
        )
        assert "MAPQ >= 25" in out
        assert "analytical 3%" in out and "calling 12%" in out
        assert "0.95" in out

    def test_says_so_when_strnaming_is_unavailable(self) -> None:
        out = render_header(RunHeader(inputs=["x.bam"], naming_markers=0))
        assert "legacy CE arithmetic" in out

    def test_no_inputs_does_not_crash(self) -> None:
        assert "none" in render_header(RunHeader())


class TestPerAlleleCoverage:
    """Integer per-allele coverage is the headline claim; it must be legible."""

    def test_genotype_carries_reads_per_allele(self) -> None:
        t = _trace(
            counts=_counts(kept=25),
            clusters=[
                _cluster(index=0, number_label="9.3", n_reads=10, called=True),
                _cluster(index=1, number_label="7", n_reads=7, called=True),
            ],
            call_rule="heterozygous",
        )
        out = render_locus(t)
        assert "9.3 (10 reads), 7 (7 reads)" in out

    def test_coverage_line_splits_called_from_discarded(self) -> None:
        t = _trace(
            counts=_counts(kept=25),
            clusters=[
                _cluster(index=0, n_reads=10, called=True),
                _cluster(index=1, n_reads=7, called=True),
                _cluster(index=2, n_reads=8, called=False),
            ],
        )
        out = render_locus(t)
        assert "25 at the locus" in out
        assert "17 on called allele(s)" in out
        assert "8 on discarded candidates" in out

    def test_haplotype_split_is_shown_per_called_allele(self) -> None:
        t = _trace(
            counts=_counts(kept=17),
            clusters=[
                _cluster(index=0, number_label="9.3", n_reads=10, n_hp1=0, n_hp2=10, called=True),
                _cluster(index=1, number_label="7", n_reads=7, n_hp1=7, n_hp2=0, called=True),
            ],
        )
        out = render_locus(t)
        assert "9.3: HP1 0 HP2 10" in out
        assert "7: HP1 7 HP2 0" in out


class TestQcVerdict:
    """The QC label is the flags themselves, never an aggregated PASS."""

    def test_a_clean_locus_says_nothing(self) -> None:
        """A green label on 95% of loci teaches a reviewer to stop reading."""
        out = render_locus(_trace(counts=_counts(), call_rule="homozygous"))
        genotype = next(x for x in out.splitlines() if "Genotype" in x)
        assert "PASS" not in genotype
        assert "OK" not in genotype

    def test_flags_are_named_on_the_genotype_line(self) -> None:
        t = _trace(counts=_counts(), flags=[("warn", "allele_imbalance")])
        genotype = next(x for x in render_locus(t).splitlines() if "Genotype" in x)
        assert "WARN" in genotype
        assert "allele_imbalance" in genotype

    def test_the_worst_severity_leads(self) -> None:
        t = _trace(
            counts=_counts(),
            flags=[("info", "hp_rescued_het"), ("warn", "low_coverage")],
        )
        genotype = next(x for x in render_locus(t).splitlines() if "Genotype" in x)
        assert genotype.index("WARN") > 0
        assert "hp_rescued_het" in genotype and "low_coverage" in genotype

    def test_allele_balance_states_its_own_scale(self) -> None:
        """0.51 alone is opaque; the reader needs to know 0.50 is even."""
        t = _trace(counts=_counts(), allele_balance=0.52)
        line = next(x for x in render_locus(t).splitlines() if "Allele balance" in x)
        assert "0.52" in line and "0.50 is even" in line and "balanced" in line

    def test_an_uneven_balance_says_uneven(self) -> None:
        t = _trace(counts=_counts(), allele_balance=0.71, balanced_ab_max=0.65)
        assert "uneven" in render_locus(t)

    def test_no_balance_line_for_a_homozygote(self) -> None:
        assert "Allele balance" not in render_locus(_trace(counts=_counts()))


class TestRunSummary:
    def test_totals_close_over_the_run(self) -> None:
        traces = [
            _trace(counts=_counts(fetched=40, kept=35), called_labels=["14"]),
            _trace(marker="vWA", counts=_counts(fetched=45, kept=38), called_labels=["14", "16"]),
        ]
        out = render_run_summary(traces)
        assert "85 fetched, 73 spanned" in out
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
