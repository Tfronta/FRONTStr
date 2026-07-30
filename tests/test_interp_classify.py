"""Tests for per-cluster classification."""

from __future__ import annotations

import pytest

from frontstr.interp.classify import classify_allele
from frontstr.interp.models import Allele, AlleleStatus


def _allele(
    *,
    consensus: str = "AGAT" * 12,
    n_reads: int = 50,
    is_deletion: bool = False,
) -> Allele:
    return Allele(
        cluster_index=0,
        consensus=consensus,
        length_bp=len(consensus),
        n_reads_total=n_reads,
        n_reads_hp1=0,
        n_reads_hp2=0,
        n_reads_hp_none=n_reads,
        n_forward=n_reads,
        n_reverse=0,
        mean_qual=30.0,
        ce=12.0,
        isfg="[AGAT]12",
        bp_diff=0,
        is_deletion=is_deletion,
    )


@pytest.mark.parametrize(
    ("n_reads", "total", "expected"),
    [
        (60, 100, AlleleStatus.ALLELE),  # 60% well above calling thresh
        (8, 100, AlleleStatus.ARTEFACT),  # 8% below calling thresh
        (1, 100, AlleleStatus.NOISE),  # 1% below analytical thresh
    ],
)
def test_classify_thresholds(n_reads: int, total: int, expected: AlleleStatus) -> None:
    a = _allele(n_reads=n_reads)
    status = classify_allele(
        a,
        total_reads=total,
        expected_stutter={},
        analytical_thresh=0.02,
        calling_thresh=0.10,
    )
    assert status == expected


def test_classify_stutter() -> None:
    """A cluster at or below its ES must be classified as stutter even if fraction is high."""
    a = _allele(consensus="AGAT" * 11, n_reads=10)
    es = {"AGAT" * 11: 12.0}
    status = classify_allele(
        a,
        total_reads=100,
        expected_stutter=es,
        analytical_thresh=0.02,
        calling_thresh=0.10,
    )
    assert status == AlleleStatus.STUTTER


def test_classify_stutter_above_es_is_allele() -> None:
    """When coverage clearly exceeds ES, status is allele."""
    a = _allele(consensus="AGAT" * 11, n_reads=50)
    es = {"AGAT" * 11: 12.0}
    status = classify_allele(
        a,
        total_reads=100,
        expected_stutter=es,
        analytical_thresh=0.02,
        calling_thresh=0.10,
    )
    assert status == AlleleStatus.ALLELE


def test_classify_deletion() -> None:
    a = _allele(consensus="", is_deletion=True, n_reads=80)
    status = classify_allele(
        a,
        total_reads=100,
        expected_stutter={},
        analytical_thresh=0.02,
        calling_thresh=0.10,
    )
    assert status == AlleleStatus.DELETION


def test_classify_no_data() -> None:
    a = _allele(n_reads=0)
    status = classify_allele(
        a,
        total_reads=0,
        expected_stutter={},
        analytical_thresh=0.02,
        calling_thresh=0.10,
    )
    assert status == AlleleStatus.NO_DATA


def test_classify_promotes_a_clean_candidate_to_allele() -> None:
    """Above both thresholds and not stutter is the only route to ALLELE.

    Used to have a sibling asserting the LongTR-driven INEXACT_ALLELE promotion;
    that path went away with LongTR (the status is retired, not repurposed).
    """
    a = _allele(n_reads=50)
    status = classify_allele(
        a,
        total_reads=100,
        expected_stutter={},
        analytical_thresh=0.02,
        calling_thresh=0.10,
    )
    assert status == AlleleStatus.ALLELE
