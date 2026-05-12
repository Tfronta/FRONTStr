"""Tests for the stutter expectation model (LUS/SLUS)."""

from __future__ import annotations

import pytest

from frontstr.evidence.cluster import Cluster
from frontstr.evidence.pileup import Observation
from frontstr.interp.stutter import (
    build_expected_stutter,
    find_motif_runs,
    top_two_runs,
    virtual_stutters,
)
from frontstr.panel.models import System


def _cluster(consensus: str, n_reads: int) -> Cluster:
    members = [
        Observation(
            read_id=f"r{i}", sequence=consensus, hp=None, mean_qual=30.0,
            strand="+", flank_left_ok=True, flank_right_ok=True,
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


def test_build_expected_stutter_default_rates() -> None:
    """100-read parent of 5 AGAT should give ES=10 for the -1 stutter (4 AGAT)."""
    system = System(
        name="TEST", chromosome="chr1", ref_start=1, ref_end=20,
        motif="AGAT", period=4,
    )
    parents = [_cluster("AGAT" * 5, n_reads=100)]
    exp = build_expected_stutter(parents, system)
    assert exp["AGAT" * 4] == pytest.approx(10.0)
    assert exp["AGAT" * 3] == pytest.approx(1.0)  # 0.10^2 * 100
    # +1 stutter uses plus_factor=0.5 by default → 5.0
    assert exp["AGAT" * 6] == pytest.approx(5.0)


def test_build_expected_stutter_overrides() -> None:
    """Per-system stutter_overrides must override defaults."""
    system = System(
        name="TPOX", chromosome="chr2", ref_start=1, ref_end=20,
        motif="AATG", period=4,
        stutter_overrides={"lus": 0.20, "plus_factor": 0.25},
    )
    parents = [_cluster("AATG" * 8, n_reads=200)]
    exp = build_expected_stutter(parents, system)
    assert exp["AATG" * 7] == pytest.approx(40.0)   # 0.20 * 200
    assert exp["AATG" * 9] == pytest.approx(10.0)   # 0.20 * 200 * 0.25


def test_build_expected_stutter_compound_motif() -> None:
    """D3-like locus: LUS-from-TCTA stutter must subtract from the longest TCTA run."""
    system = System(
        name="D3S1358", chromosome="chr3", ref_start=1, ref_end=100,
        motif="TCTA,TCTG", period=4,
    )
    parent_seq = "TCTA" * 5 + "TCTG" + "TCTA" * 8
    parents = [_cluster(parent_seq, n_reads=100)]
    exp = build_expected_stutter(parents, system)
    # The -1 LUS variant collapses the 8-TCTA run to 7
    minus_one_lus = "TCTA" * 5 + "TCTG" + "TCTA" * 7
    assert minus_one_lus in exp
    assert exp[minus_one_lus] == pytest.approx(10.0)
    # The -1 SLUS variant collapses the 5-TCTA run to 4 (slus_rate=0.05)
    minus_one_slus = "TCTA" * 4 + "TCTG" + "TCTA" * 8
    assert minus_one_slus in exp
    assert exp[minus_one_slus] == pytest.approx(5.0)


def test_build_expected_stutter_accumulates() -> None:
    """Two parents that share a -1 variant should accumulate their contributions."""
    system = System(
        name="T", chromosome="chr1", ref_start=1, ref_end=100,
        motif="AGAT", period=4,
    )
    parents = [
        _cluster("AGAT" * 5, n_reads=100),
        _cluster("AGAT" * 5, n_reads=50),
    ]
    exp = build_expected_stutter(parents, system)
    assert exp["AGAT" * 4] == pytest.approx(15.0)


def test_build_expected_stutter_empty() -> None:
    system = System(
        name="T", chromosome="chr1", ref_start=1, ref_end=10,
        motif="AGAT", period=4,
    )
    assert build_expected_stutter([], system) == {}
