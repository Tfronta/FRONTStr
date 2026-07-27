"""Tests for the stutter expectation model (LUS/SLUS)."""

from __future__ import annotations

import math

import pytest

from frontstr.evidence.cluster import Cluster
from frontstr.evidence.pileup import Observation
from frontstr.interp.stutter import (
    DEFAULT_STUTTER_MODEL,
    build_expected_stutter,
    find_motif_runs,
    top_two_runs,
    virtual_stutters,
)
from frontstr.panel.models import System


def _cluster(consensus: str, n_reads: int) -> Cluster:
    members = [
        Observation(
            read_id=f"r{i}",
            sequence=consensus,
            hp=None,
            mean_qual=30.0,
            strand="+",
            flank_left_ok=True,
            flank_right_ok=True,
        )
        for i in range(n_reads)
    ]
    return Cluster(consensus=consensus, members=members)


def test_find_motif_runs_simple() -> None:
    """A simple repeat: 5 copies of AGAT."""
    seq = "AGAT" * 5
    runs = find_motif_runs(seq, ["AGAT"])
    assert len(runs) == 1
    assert runs[0].motif == "AGAT"
    assert runs[0].start == 0
    assert runs[0].n_copies == 5


def test_find_motif_runs_compound() -> None:
    """D3S1358-like: 5 TCTA, 1 TCTG, 4 TCTA."""
    seq = "TCTA" * 5 + "TCTG" + "TCTA" * 4
    runs = find_motif_runs(seq, ["TCTA", "TCTG"])
    # Greedy scan: first 5 TCTA, then 1 TCTG, then 4 TCTA
    assert [r.n_copies for r in runs] == [5, 1, 4]
    assert [r.motif for r in runs] == ["TCTA", "TCTG", "TCTA"]


def test_top_two_runs() -> None:
    seq = "TCTA" * 6 + "TCTG" + "TCTA" * 3
    runs = find_motif_runs(seq, ["TCTA", "TCTG"])
    lus, slus = top_two_runs(runs)
    assert lus is not None and lus.n_copies == 6
    assert slus is not None and slus.n_copies == 3


def test_top_two_runs_empty() -> None:
    assert top_two_runs([]) == (None, None)


def test_virtual_stutters_minus_one() -> None:
    """Parent with 5 AGAT must yield a -1 LUS variant of 4 AGAT."""
    seq = "AGAT" * 5
    variants = list(virtual_stutters(seq, ["AGAT"]))
    targets = {v[0] for v in variants}
    assert ("AGAT" * 4) in targets
    assert ("AGAT" * 3) in targets  # -2
    assert ("AGAT" * 6) in targets  # +1


def test_virtual_stutters_skip_singleton_runs() -> None:
    """A 1-copy run cannot produce a stutter without disappearing."""
    seq = "TCTA"
    variants = list(virtual_stutters(seq, ["TCTA"]))
    assert variants == []


def test_build_expected_stutter_scales_with_lus() -> None:
    """The rate is a function of the slipping run's length, not a constant.

    This is the whole point of the calibrated model: a flat per-marker rate
    cannot express that stutter at LUS 14 is ~6x that at LUS 10 (measured on
    the ONT R10 slice set).
    """
    system = System(
        name="TEST",
        chromosome="chr1",
        ref_start=1,
        ref_end=200,
        motif="AGAT",
        period=4,
    )
    short = build_expected_stutter([_cluster("AGAT" * 10, n_reads=100)], system)
    long = build_expected_stutter([_cluster("AGAT" * 14, n_reads=100)], system)
    assert long["AGAT" * 13] > 4 * short["AGAT" * 9]


def test_build_expected_stutter_step_factors() -> None:
    """-2 and +1 are multipliers of the -1 rate, not a geometric decay.

    The geometric ``rate ** step`` form underestimates -2 on ONT by ~2.5x.
    """
    system = System(
        name="TEST",
        chromosome="chr1",
        ref_start=1,
        ref_end=200,
        motif="AGAT",
        period=4,
    )
    exp = build_expected_stutter([_cluster("AGAT" * 13, n_reads=100)], system)
    minus1 = exp["AGAT" * 12]
    factors = DEFAULT_STUTTER_MODEL.step_factors
    assert exp["AGAT" * 11] == pytest.approx(minus1 * factors["-2"])
    assert exp["AGAT" * 14] == pytest.approx(minus1 * factors["1"])
    # -2 is a multiplier of -1, NOT the geometric rate**2 the old model used.
    assert factors["-2"] > minus1 / 100 * 2


def test_lus_is_clamped_to_the_calibrated_range() -> None:
    """Outside the fitted LUS range, use the nearest supported rate.

    Extrapolating the line would predict a negative (i.e. zero) rate at low LUS
    and an implausible one at very high LUS; neither is measured.
    """
    system = System(
        name="TEST",
        chromosome="chr1",
        ref_start=1,
        ref_end=400,
        motif="AGAT",
        period=4,
    )
    at_min = build_expected_stutter([_cluster("AGAT" * 10, n_reads=100)], system)
    below = build_expected_stutter([_cluster("AGAT" * 6, n_reads=100)], system)
    assert below["AGAT" * 5] == pytest.approx(at_min["AGAT" * 9])

    at_max = build_expected_stutter([_cluster("AGAT" * 15, n_reads=100)], system)
    above = build_expected_stutter([_cluster("AGAT" * 30, n_reads=100)], system)
    assert above["AGAT" * 29] == pytest.approx(at_max["AGAT" * 14])


def test_slus_stutter_is_rated_by_its_own_run_length() -> None:
    """A shorter secondary run must produce less stutter than the primary one."""
    system = System(
        name="D3S1358",
        chromosome="chr3",
        ref_start=1,
        ref_end=200,
        motif="TCTA,TCTG",
        period=4,
    )
    parent_seq = "TCTA" * 11 + "TCTG" + "TCTA" * 14
    exp = build_expected_stutter([_cluster(parent_seq, n_reads=100)], system)

    minus_one_lus = "TCTA" * 11 + "TCTG" + "TCTA" * 13
    minus_one_slus = "TCTA" * 10 + "TCTG" + "TCTA" * 14
    assert exp[minus_one_lus] > exp[minus_one_slus] > 0


def test_legacy_flat_lus_override_is_still_honoured() -> None:
    """Labs validated against the old flat rate can pin it per marker."""
    system = System(
        name="TPOX",
        chromosome="chr2",
        ref_start=1,
        ref_end=200,
        motif="AATG",
        period=4,
        stutter_overrides={"lus": 0.20},
    )
    exp = build_expected_stutter([_cluster("AATG" * 8, n_reads=200)], system)
    assert exp["AATG" * 7] == pytest.approx(40.0)  # 0.20 * 200


def test_per_marker_slope_override() -> None:
    system = System(
        name="TEST",
        chromosome="chr1",
        ref_start=1,
        ref_end=200,
        motif="AGAT",
        period=4,
        stutter_overrides={"log_intercept": -4.0, "log_slope": 0.0},
    )
    exp = build_expected_stutter([_cluster("AGAT" * 12, n_reads=100)], system)
    assert exp["AGAT" * 11] == pytest.approx(math.exp(-4.0) * 100)


def test_build_expected_stutter_accumulates() -> None:
    """Two parents that share a -1 variant should accumulate their contributions."""
    system = System(
        name="T",
        chromosome="chr1",
        ref_start=1,
        ref_end=200,
        motif="AGAT",
        period=4,
    )
    one = build_expected_stutter([_cluster("AGAT" * 12, n_reads=100)], system)
    two = build_expected_stutter(
        [_cluster("AGAT" * 12, n_reads=100), _cluster("AGAT" * 12, n_reads=50)], system
    )
    assert two["AGAT" * 11] == pytest.approx(one["AGAT" * 11] * 1.5)


def test_build_expected_stutter_empty() -> None:
    system = System(
        name="T",
        chromosome="chr1",
        ref_start=1,
        ref_end=10,
        motif="AGAT",
        period=4,
    )
    assert build_expected_stutter([], system) == {}
