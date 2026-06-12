"""Tests for tri-allelic / mixture call profile."""

from __future__ import annotations

from frontstr.interp.models import Allele, AlleleStatus, CallRule, TriType
from frontstr.interp.triallelic import call_profile
from frontstr.panel.models import System


def _allele(idx: int, n_reads: int, *, status: AlleleStatus = AlleleStatus.ALLELE,
            consensus: str | None = None) -> Allele:
    return Allele(
        cluster_index=idx,
        consensus=consensus or f"SEQ{idx}",
        length_bp=10 + idx,
        n_reads_total=n_reads, n_reads_hp1=0, n_reads_hp2=0, n_reads_hp_none=n_reads,
        n_forward=n_reads, n_reverse=0, mean_qual=30.0,
        ce=float(10 + idx), isfg="[X]10", bp_diff=0,
        is_deletion=False, status=status,
    )


def _system(allow_tri: bool = False, tri_balanced_thr: float | None = None,
            name: str = "TPOX") -> System:
    return System(
        name=name, chromosome="chr2", ref_start=1_490_000, ref_end=1_490_040,
        motif="AATG", period=4,
        allow_triallelic=allow_tri, tri_balanced_thr=tri_balanced_thr,
    )


def test_no_candidates_no_data() -> None:
    sys = _system()
    a, rule, tri = call_profile([], sys)
    assert a == [] and rule == CallRule.NO_DATA and tri == TriType.NONE


def test_single_allele_homozygous() -> None:
    sys = _system()
    a, rule, _tri = call_profile([_allele(0, 80)], sys)
    assert len(a) == 1 and rule == CallRule.HOMOZYGOUS


def test_low_second_allele_collapses_to_hom() -> None:
    sys = _system()
    a1 = _allele(0, 100)
    a2 = _allele(1, 25)  # 25% - below default 40% het floor
    a, rule, _tri = call_profile([a1, a2], sys)
    assert a == [a1] and rule == CallRule.HOMOZYGOUS


def test_balanced_pair_heterozygous() -> None:
    sys = _system()
    a, rule, tri = call_profile([_allele(0, 60), _allele(1, 55)], sys)
    assert len(a) == 2 and rule == CallRule.HETEROZYGOUS and tri == TriType.NONE


def test_three_peaks_in_non_triallelic_locus_flags_mixture() -> None:
    """vWA-like locus without allow_triallelic — third peak triggers mixture review."""
    sys = _system(allow_tri=False)
    a, rule, tri = call_profile(
        [_allele(0, 60), _allele(1, 55), _allele(2, 50)], sys
    )
    assert len(a) == 2  # only top two reported
    assert rule == CallRule.TWO_CALLED_THREE_PRESENT
    assert tri == TriType.MIXTURE_SUSPECTED


def test_tpox_type_ii_balanced() -> None:
    """TPOX with allow_triallelic + tri_balanced_thr=0.5: three balanced peaks → Type II."""
    sys = _system(allow_tri=True, tri_balanced_thr=0.5)
    a, rule, tri = call_profile(
        [_allele(0, 60), _allele(1, 50), _allele(2, 45)], sys
    )
    assert len(a) == 3
    assert rule == CallRule.TRIALLELIC_TYPE_II
    assert tri == TriType.TYPE_II_BALANCED


def test_tpox_type_i_unbalanced() -> None:
    """Top two balanced + small third → Type I."""
    sys = _system(allow_tri=True, tri_balanced_thr=0.6)
    a, rule, tri = call_profile(
        [_allele(0, 100), _allele(1, 90), _allele(2, 35)], sys
    )
    assert len(a) == 3
    assert rule == CallRule.TRIALLELIC_TYPE_I
    assert tri == TriType.TYPE_I_UNBALANCED


def test_triallelic_review_when_pattern_is_weird() -> None:
    """Top two unbalanced but 3rd present — needs review."""
    sys = _system(allow_tri=True, tri_balanced_thr=0.6)
    # phr_12 = 0.4 (below balanced threshold)
    _a, rule, tri = call_profile(
        [_allele(0, 100), _allele(1, 40), _allele(2, 35)], sys
    )
    assert rule == CallRule.TRIALLELIC_REVIEW
    assert tri == TriType.MIXTURE_SUSPECTED


def test_third_below_calling_thresh_drops_to_het() -> None:
    """A 3rd allele below calling threshold should not promote to tri."""
    sys = _system(allow_tri=True)
    # Third is 5% of top - below calling_thresh of 10% by default
    a, rule, _tri = call_profile(
        [_allele(0, 100), _allele(1, 80), _allele(2, 5)], sys,
        calling_thresh=0.10,
    )
    assert len(a) == 2 and rule == CallRule.HETEROZYGOUS


def test_third_below_read_floor_drops_to_het() -> None:
    """Bug #6: a 3rd candidate above calling_thresh fraction but below the
    absolute read floor (default 5) is an ONT phantom — call het, no mixture."""
    sys = _system(allow_tri=False)  # default min_reads_third = 5
    a1, a2, a3 = _allele(0, 30), _allele(1, 28), _allele(2, 4)
    # phr_13 = 4/30 = 0.133 >= 0.10 (passes the fractional gate) but 4 < 5 floor
    a, rule, tri = call_profile([a1, a2, a3], sys, calling_thresh=0.10)
    assert a == [a1, a2]
    assert rule == CallRule.HETEROZYGOUS
    assert tri == TriType.NONE


def test_read_floor_configurable_via_system() -> None:
    """A higher per-marker min_reads_third floor suppresses larger phantoms."""
    sys = System(
        name="D8S1179", chromosome="chr8", ref_start=1, ref_end=40,
        motif="TCTA", period=4, allow_triallelic=False, min_reads_third=10,
    )
    a1, a2, a3 = _allele(0, 60), _allele(1, 55), _allele(2, 8)
    # 8/60 = 0.133 passes fraction, but 8 < 10 floor → het, not mixture
    a, rule, tri = call_profile([a1, a2, a3], sys, calling_thresh=0.10)
    assert a == [a1, a2]
    assert rule == CallRule.HETEROZYGOUS
    assert tri == TriType.NONE


def test_real_third_above_floor_still_flags_mixture() -> None:
    """A 3rd allele clearing both the fraction and the read floor still flags."""
    sys = _system(allow_tri=False)
    a, rule, tri = call_profile(
        [_allele(0, 60), _allele(1, 55), _allele(2, 50)], sys
    )
    assert len(a) == 2
    assert rule == CallRule.TWO_CALLED_THREE_PRESENT
    assert tri == TriType.MIXTURE_SUSPECTED


def test_only_allele_status_considered() -> None:
    """Stutter/noise candidates are excluded before counting."""
    sys = _system(allow_tri=True)
    a, rule, _tri = call_profile(
        [
            _allele(0, 100, status=AlleleStatus.ALLELE),
            _allele(1, 80, status=AlleleStatus.STUTTER),
            _allele(2, 70, status=AlleleStatus.ALLELE),
        ],
        sys,
    )
    # only the two ALLELEs survive
    assert len(a) == 2 and rule == CallRule.HETEROZYGOUS
