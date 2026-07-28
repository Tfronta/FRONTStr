"""Tests for :mod:`frontstr.interp.haplotype` — same-haplotype phantom suppression.

The contract: a candidate is suppressed only when the haplotype tags prove it
is the same allele as a stronger cluster. Everything else — unphased data,
weak tagging, a whole repeat unit of separation, a different haplotype, a
triallelic-capable marker — must leave the candidate standing.
"""

from __future__ import annotations

import pytest

from frontstr.interp.haplotype import (
    DEFAULT_MAX_PHANTOM_BP,
    dominant_hp,
    on_opposite_haplotypes,
    suppress_hp_phantoms,
)
from frontstr.interp.models import Allele, AlleleStatus
from frontstr.panel.models import System


def _system(*, allow_triallelic: bool = False) -> System:
    return System(
        name="D8S1179",
        chromosome="chr8",
        ref_start=1,
        ref_end=300,
        motif="TCTA",
        period=4,
        allow_triallelic=allow_triallelic,
    )


def _allele(
    idx: int,
    length_bp: int,
    *,
    hp1: int = 0,
    hp2: int = 0,
    untagged: int = 0,
    status: AlleleStatus = AlleleStatus.ALLELE,
) -> Allele:
    n = hp1 + hp2 + untagged
    return Allele(
        cluster_index=idx,
        consensus="A" * length_bp,
        length_bp=length_bp,
        n_reads_total=n,
        n_reads_hp1=hp1,
        n_reads_hp2=hp2,
        n_reads_hp_none=untagged,
        n_forward=n // 2,
        n_reverse=n - n // 2,
        mean_qual=30.0,
        ce=None,
        isfg="",
        bp_diff=0,
        is_deletion=False,
        status=status,
    )


def _statuses(alleles: list[Allele]) -> list[str]:
    return [a.status.value for a in alleles]


# ---------------------------------------------------------------------------
# dominant_hp
# ---------------------------------------------------------------------------


def test_untagged_cluster_has_no_haplotype() -> None:
    assert dominant_hp(_allele(0, 250, untagged=12)) is None


def test_too_few_tagged_reads_is_not_evidence() -> None:
    """Two tagged reads can agree by chance; that is not a haplotype claim."""
    assert dominant_hp(_allele(0, 250, hp1=2, untagged=8)) is None
    assert dominant_hp(_allele(0, 250, hp1=3, untagged=8)) == 1


def test_untagged_reads_do_not_dilute_purity() -> None:
    """whatshap leaving a read unplaced is missing evidence, not counter-evidence."""
    assert dominant_hp(_allele(0, 250, hp1=5, untagged=20)) == 1


def test_mixed_haplotype_cluster_is_unassigned() -> None:
    assert dominant_hp(_allele(0, 250, hp1=5, hp2=4)) is None


# ---------------------------------------------------------------------------
# Suppression fires
# ---------------------------------------------------------------------------


def test_same_length_same_haplotype_is_collapsed() -> None:
    """The identity-split signature: same bp, same HP, split by the 0.97 threshold."""
    owner = _allele(0, 252, hp1=13)
    phantom = _allele(1, 252, hp1=5)
    other = _allele(2, 240, hp2=12)

    assert suppress_hp_phantoms([owner, phantom, other], _system()) == 1
    assert phantom.status == AlleleStatus.HP_PHANTOM
    assert owner.status == AlleleStatus.ALLELE
    assert other.status == AlleleStatus.ALLELE


def test_one_bp_apart_same_haplotype_is_collapsed() -> None:
    """The length-split signature: a 1 bp indel error crossing the length bin."""
    owner = _allele(0, 289, hp1=10)
    phantom = _allele(1, 288, hp1=6)
    assert suppress_hp_phantoms([owner, phantom], _system()) == 1
    assert phantom.status == AlleleStatus.HP_PHANTOM


def test_absorbed_reads_are_recorded_on_the_owner() -> None:
    """Coverage must not be understated once a phantom is folded in."""
    owner = _allele(0, 252, hp1=13)
    phantom = _allele(1, 252, hp1=5)
    suppress_hp_phantoms([owner, phantom], _system())
    assert owner.n_reads_absorbed == 5
    assert owner.n_reads_total == 13, "observed count stays as observed"


def test_strongest_cluster_wins_ownership() -> None:
    weak = _allele(0, 252, hp1=4)
    strong = _allele(1, 252, hp1=11)
    suppress_hp_phantoms([weak, strong], _system())
    assert strong.status == AlleleStatus.ALLELE
    assert weak.status == AlleleStatus.HP_PHANTOM
    assert strong.n_reads_absorbed == 4


# ---------------------------------------------------------------------------
# Suppression must NOT fire
# ---------------------------------------------------------------------------


def test_unphased_data_is_a_no_op() -> None:
    alleles = [_allele(0, 252, untagged=14), _allele(1, 252, untagged=5)]
    assert suppress_hp_phantoms(alleles, _system()) == 0
    assert _statuses(alleles) == ["allele", "allele"]


def test_different_haplotypes_are_never_collapsed() -> None:
    """Same length on opposite haplotypes is a real homozygote-by-length, not a phantom."""
    alleles = [_allele(0, 252, hp1=12), _allele(1, 252, hp2=10)]
    assert suppress_hp_phantoms(alleles, _system()) == 0
    assert _statuses(alleles) == ["allele", "allele"]


def test_a_full_repeat_unit_apart_is_a_different_allele() -> None:
    """4 bp on a tetramer is a real neighbouring allele — never absorb it."""
    alleles = [_allele(0, 252, hp1=12), _allele(1, 248, hp1=6)]
    assert suppress_hp_phantoms(alleles, _system()) == 0
    assert _statuses(alleles) == ["allele", "allele"]


@pytest.mark.parametrize("delta", [DEFAULT_MAX_PHANTOM_BP, DEFAULT_MAX_PHANTOM_BP + 1])
def test_max_phantom_bp_boundary(delta: int) -> None:
    alleles = [_allele(0, 252, hp1=12), _allele(1, 252 - delta, hp1=6)]
    collapsed = suppress_hp_phantoms(alleles, _system())
    assert collapsed == (1 if delta <= DEFAULT_MAX_PHANTOM_BP else 0)


def test_triallelic_marker_is_skipped() -> None:
    """A duplication genuinely puts two alleles on one haplotype."""
    alleles = [_allele(0, 252, hp1=12), _allele(1, 252, hp1=6)]
    assert suppress_hp_phantoms(alleles, _system(allow_triallelic=True)) == 0
    assert _statuses(alleles) == ["allele", "allele"]


def test_non_candidate_statuses_are_ignored() -> None:
    """Stutter and noise clusters were already resolved; do not re-litigate them."""
    owner = _allele(0, 252, hp1=12)
    stutter = _allele(1, 252, hp1=5, status=AlleleStatus.STUTTER)
    assert suppress_hp_phantoms([owner, stutter], _system()) == 0
    assert stutter.status == AlleleStatus.STUTTER
    assert owner.n_reads_absorbed == 0


# ---------------------------------------------------------------------------
# The case a read-count floor gets wrong (known-bug #6, GM19038 D12S391)
# ---------------------------------------------------------------------------


def test_weak_real_allele_survives_when_the_strong_pair_are_both_phantoms() -> None:
    """Two HP1 clusters plus a weaker genuine HP2 allele.

    A read-count floor would drop the 5-read HP2 cluster and emit a confidently
    wrong homozygote. Haplotype evidence keeps it and removes the right one.
    """
    hp1_owner = _allele(0, 289, hp1=10)
    hp1_phantom = _allele(1, 288, hp1=6)
    hp2_real = _allele(2, 301, hp2=5)
    hp2_phantom = _allele(3, 302, hp2=4)

    alleles = [hp1_owner, hp1_phantom, hp2_real, hp2_phantom]
    assert suppress_hp_phantoms(alleles, _system()) == 2

    survivors = [a for a in alleles if a.status == AlleleStatus.ALLELE]
    assert [a.length_bp for a in survivors] == [289, 301]
    assert hp1_owner.n_reads_absorbed == 6
    assert hp2_real.n_reads_absorbed == 4


# ---------------------------------------------------------------------------
# on_opposite_haplotypes — the same invariant, read backwards
# ---------------------------------------------------------------------------


def test_opposite_haplotypes_detected() -> None:
    assert on_opposite_haplotypes(_allele(0, 100, hp1=17), _allele(1, 104, hp2=5)) is True


def test_same_haplotype_is_not_opposite() -> None:
    assert on_opposite_haplotypes(_allele(0, 100, hp1=17), _allele(1, 104, hp2=0, hp1=5)) is False


def test_unphased_pair_is_not_opposite() -> None:
    """No tags at all — the whole mechanism must be a no-op."""
    assert (
        on_opposite_haplotypes(_allele(0, 100, untagged=17), _allele(1, 104, untagged=5)) is False
    )


def test_one_side_below_the_tagged_floor_is_not_opposite() -> None:
    assert on_opposite_haplotypes(_allele(0, 100, hp1=17), _allele(1, 104, hp2=2)) is False


def test_impure_side_is_not_opposite() -> None:
    """A cluster straddling both haplotypes cannot anchor the invariant."""
    assert on_opposite_haplotypes(_allele(0, 100, hp1=11, hp2=7), _allele(1, 104, hp2=5)) is False


def test_untagged_reads_do_not_block_the_call() -> None:
    """Missing phasing evidence is not evidence against, as in dominant_hp."""
    a = _allele(0, 100, hp1=9, untagged=8)
    b = _allele(1, 104, hp2=4, untagged=3)
    assert on_opposite_haplotypes(a, b) is True
