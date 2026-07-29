"""Tests for the run-parameter table.

The table exists so a run can say what it used and so a laboratory can change a
threshold to test one. These tests guard the part that makes both safe: that
overriding a *measured* default is distinguishable from ordinary tuning.
"""

from __future__ import annotations

import pytest

from frontstr.params import PARAM_SPECS, RunParameters, render_echo


def test_every_spec_has_a_reason() -> None:
    """A default nobody can justify is a magic number with better manners."""
    for spec in PARAM_SPECS:
        assert spec.note.strip(), f"{spec.name} has no note"
        assert spec.provenance in ("derived", "chosen", "convention")


def test_names_are_unique() -> None:
    names = [s.name for s in PARAM_SPECS]
    assert len(names) == len(set(names))


def test_defaults_match_the_code_they_describe() -> None:
    """The table is documentation until it agrees with the implementation.

    If one of these drifts, a report will confidently print a default the
    pipeline is not using.
    """
    from frontstr.evidence.cluster import _DEFAULT_IDENTITY_THRESHOLD
    from frontstr.evidence.pileup import _DEFAULT_FLANK_ANCHOR
    from frontstr.interp.profile import DEFAULT_ANALYTICAL_THRESH, DEFAULT_CALLING_THRESH
    from frontstr.interp.qc import QcThresholds
    from frontstr.interp.triallelic import DEFAULT_MIN_PHR_FOR_HET
    from frontstr.panel.models import System

    params = RunParameters.defaults()
    qc = QcThresholds()
    assert params["flank_anchor"] == _DEFAULT_FLANK_ANCHOR
    assert params["identity_threshold"] == _DEFAULT_IDENTITY_THRESHOLD
    assert params["analytical_thresh"] == DEFAULT_ANALYTICAL_THRESH
    assert params["calling_thresh"] == DEFAULT_CALLING_THRESH
    assert params["min_phr_for_het"] == DEFAULT_MIN_PHR_FOR_HET
    assert params["low_coverage_reads"] == qc.low_coverage_reads
    assert params["balanced_ab_max"] == qc.balanced_ab_max
    assert params["min_reads_third"] == System.model_fields["min_reads_third"].default


def test_of_fills_unset_values_with_defaults() -> None:
    p = RunParameters.of(min_mapq=30)
    assert p["min_mapq"] == 30
    assert p["calling_thresh"] == 0.10
    assert not p.is_default("min_mapq")
    assert p.is_default("calling_thresh")


def test_none_means_unset_not_null() -> None:
    """Typer passes None for an unsupplied option; that must mean 'default'."""
    assert RunParameters.of(min_reads_third=None)["min_reads_third"] == 5


def test_unknown_parameter_is_rejected() -> None:
    """A typo must fail loudly rather than be silently ignored."""
    with pytest.raises(KeyError, match="min_maqp"):
        RunParameters.of(min_maqp=30)


class TestOverrides:
    def test_a_default_run_has_none(self) -> None:
        assert RunParameters.defaults().overrides() == []
        assert RunParameters.defaults().derived_overrides() == []

    def test_tuning_a_chosen_threshold_does_not_mark_the_run(self) -> None:
        """Changing calling_thresh is ordinary work, not a provenance event."""
        p = RunParameters.of(calling_thresh=0.05)
        assert [s.name for s in p.overrides()] == ["calling_thresh"]
        assert p.derived_overrides() == []

    def test_lowering_a_measured_default_marks_the_run(self) -> None:
        p = RunParameters.of(min_reads_third=2)
        assert [s.name for s in p.derived_overrides()] == ["min_reads_third"]

    def test_setting_a_default_to_its_own_value_is_not_an_override(self) -> None:
        assert RunParameters.of(min_reads_third=5).overrides() == []


class TestEcho:
    def test_shows_every_parameter_not_only_the_changed_ones(self) -> None:
        """A run listing only what was typed reads as barely configured."""
        out = render_echo(RunParameters.defaults())
        for spec in PARAM_SPECS:
            assert spec.name in out

    def test_marks_a_changed_value_with_its_default(self) -> None:
        out = render_echo(RunParameters.of(min_mapq=5))
        line = next(x for x in out.splitlines() if "min_mapq" in x)
        assert "CHANGED" in line
        assert "default 20" in line

    def test_explains_why_a_measured_default_mattered(self) -> None:
        out = render_echo(RunParameters.of(min_reads_third=1))
        assert "not comparable with a default one" in out
        assert "known-bug #6" in out, "the reason must travel with the warning"

    def test_a_default_run_says_nothing_alarming(self) -> None:
        out = render_echo(RunParameters.defaults())
        assert "CHANGED" not in out
        assert "not comparable" not in out


class TestAuditRows:
    def test_one_row_per_parameter_with_provenance(self) -> None:
        rows = RunParameters.of(min_reads_third=2).as_audit_rows()
        assert len(rows) == len(PARAM_SPECS)
        row = next(r for r in rows if r["name"] == "min_reads_third")
        assert row["value"] == 2
        assert row["default"] == 5
        assert row["is_default"] is False
        assert row["provenance"] == "derived"
        assert row["note"]

    def test_unchanged_rows_are_marked_default(self) -> None:
        rows = RunParameters.defaults().as_audit_rows()
        assert all(r["is_default"] for r in rows)
