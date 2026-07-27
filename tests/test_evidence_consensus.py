"""Tests for :mod:`frontstr.evidence.consensus` — the POA backend chain.

The forensic contract under test: a cluster's consensus must be a *polished*
multiple-alignment consensus whenever a POA backend is installed, and the
record must always say which method produced it so an unpolished consensus is
never mistaken for a polished one.
"""

from __future__ import annotations

import random

import pytest

from frontstr.evidence.consensus import (
    ConsensusMethod,
    build_consensus,
    poa_backend_name,
    reset_backend_cache,
)

poa_only = pytest.mark.skipif(
    not poa_backend_name(), reason="no POA backend installed (pyabpoa / pyspoa)"
)


@pytest.fixture(autouse=True)
def _restore_backend_cache():
    """Every test starts from a clean backend resolution."""
    reset_backend_cache()
    yield
    reset_backend_cache()


@pytest.fixture
def no_poa(monkeypatch: pytest.MonkeyPatch) -> None:
    """Force the mode fallback, reproducing an install with no POA backend."""
    monkeypatch.setattr("frontstr.evidence.consensus._BACKEND", None)


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------


def test_empty_cluster_is_reported_as_empty() -> None:
    assert build_consensus([]) == ("", ConsensusMethod.EMPTY)


def test_single_read_is_returned_verbatim() -> None:
    """One read cannot be polished; the method must say so rather than lie."""
    seq, method = build_consensus(["AGATAGATAGAT"])
    assert seq == "AGATAGATAGAT"
    assert method == ConsensusMethod.SINGLE
    assert not method.is_poa


def test_identical_reads_agree_under_every_backend() -> None:
    seq, method = build_consensus(["AGAT" * 6] * 5)
    assert seq == "AGAT" * 6
    assert method in (
        ConsensusMethod.POA_ABPOA,
        ConsensusMethod.POA_SPOA,
        ConsensusMethod.MODE,
    )


# ---------------------------------------------------------------------------
# Fallback behaviour
# ---------------------------------------------------------------------------


def test_without_backend_falls_back_to_mode(no_poa: None) -> None:
    seqs = ["AGATAGAT", "AGATAGAT", "AGATAGAT", "AGATTGAT"]
    seq, method = build_consensus(seqs)
    assert seq == "AGATAGAT"
    assert method == ConsensusMethod.MODE
    assert not method.is_poa


def test_backend_failure_degrades_instead_of_losing_the_locus(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A POA crash must not drop the marker — it degrades to mode, flagged."""

    class _Exploding:
        method = ConsensusMethod.POA_SPOA

        def consensus(self, seqs: list[str]) -> str:
            raise RuntimeError("boom")

    monkeypatch.setattr("frontstr.evidence.consensus._BACKEND", _Exploding())
    seq, method = build_consensus(["AGATAGAT"] * 3 + ["AGATTGAT"])
    assert seq == "AGATAGAT"
    assert method == ConsensusMethod.MODE


def test_backend_returning_empty_degrades_to_mode(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _Silent:
        method = ConsensusMethod.POA_SPOA

        def consensus(self, seqs: list[str]) -> str:
            return ""

    monkeypatch.setattr("frontstr.evidence.consensus._BACKEND", _Silent())
    seq, method = build_consensus(["AGATAGAT"] * 3)
    assert seq == "AGATAGAT"
    assert method == ConsensusMethod.MODE


def test_backend_name_is_empty_without_a_backend(no_poa: None) -> None:
    assert poa_backend_name() == ""


# ---------------------------------------------------------------------------
# The reason this module exists: POA corrects errors that mode cannot
# ---------------------------------------------------------------------------


def _mutate(seq: str, rng: random.Random, sub: float = 0.01) -> str:
    """Apply substitution-only ONT-like noise, preserving length.

    Length is preserved so the clustering contract (length-binned members) is
    respected and the comparison isolates base-level correction.
    """
    return "".join(
        rng.choice([b for b in "ACGT" if b != c]) if rng.random() < sub else c
        for c in seq
    )


@poa_only
def test_poa_recovers_the_true_haplotype_where_mode_cannot() -> None:
    """The measured claim from the module docstring, as a regression test.

    With independent per-read errors, the most common *exact* sequence is only
    the truth when some read happens to be error-free — a coin flip that gets
    worse as the locus gets longer. POA votes per column and recovers the truth
    regardless. This is what makes iso-allele calling trustworthy.
    """
    rng = random.Random(11)
    truth = (
        "".join(rng.choice("ACGT") for _ in range(100))
        + "TCTA" * 13
        + "".join(rng.choice("ACGT") for _ in range(100))
    )
    reads = [_mutate(truth, rng, sub=0.02) for _ in range(12)]
    # Precondition: no read is error-free, so the mode *cannot* be the truth.
    assert truth not in reads

    poa_seq, method = build_consensus(reads)
    assert method.is_poa
    assert poa_seq == truth

    # Same reads, no backend: the mode is one of the noisy reads, not the truth.
    from frontstr.evidence import consensus as mod

    mod._BACKEND = None
    mode_seq, mode_method = build_consensus(reads)
    assert mode_method == ConsensusMethod.MODE
    assert mode_seq != truth


@poa_only
def test_poa_preserves_length_of_a_length_binned_cluster() -> None:
    """Global alignment must not trim flanks — length *is* the CE number.

    Local alignment would be free to clip the flanks, silently changing the
    called allele length and therefore its CE designation.
    """
    rng = random.Random(3)
    truth = "ACGTACGTGGCC" + "AGAT" * 11 + "TTGCAATTGCAA"
    reads = [_mutate(truth, rng, sub=0.02) for _ in range(15)]
    seq, method = build_consensus(reads)
    assert method.is_poa
    assert len(seq) == len(truth)


@poa_only
def test_poa_does_not_invent_a_repeat_unit_in_the_flank() -> None:
    """Read errors that spell the motif in a flank must not survive polishing.

    Observed on real data (HG00263 D18S51): the mode consensus contained a
    spurious ``AGAA`` in the left flank, which shifted the bracket structure.
    """
    rng = random.Random(5)
    truth = "GAGGCAGGAGGAGTTCTTGAGCCC" + "AGAA" * 12 + "GTTAATTTTAATTTTAACATGT"
    reads = [_mutate(truth, rng, sub=0.015) for _ in range(14)]
    seq, method = build_consensus(reads)
    assert method.is_poa
    assert seq == truth
