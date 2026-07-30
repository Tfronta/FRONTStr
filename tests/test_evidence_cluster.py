"""Tests for :func:`frontstr.evidence.cluster.cluster_observations`."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.errors import EvidenceError
from frontstr.evidence.cluster import Cluster, cluster_observations
from frontstr.evidence.pileup import Observation, pileup_locus
from tests.conftest import SYNTH_CHROM, SYNTH_TR_END, SYNTH_TR_START


def _obs(
    sequence: str, hp: int | None = None, name: str = "r", ps: int | None = None
) -> Observation:
    return Observation(
        read_id=name,
        sequence=sequence,
        hp=hp,
        ps=ps,
        mean_qual=40.0,
        strand="+",
        flank_left_ok=True,
        flank_right_ok=True,
    )


def test_empty_input() -> None:
    assert cluster_observations([]) == []


def test_singleton_cluster() -> None:
    out = cluster_observations([_obs("AGAT" * 5)])
    assert len(out) == 1
    assert out[0].n_reads == 1
    assert out[0].consensus == "AGAT" * 5


def test_identical_sequences_collapse() -> None:
    obs = [_obs("AGAT" * 5, hp=1, name=f"r{i}") for i in range(5)]
    out = cluster_observations(obs)
    assert len(out) == 1
    assert out[0].n_reads == 5
    assert out[0].n_hp1 == 5


def test_different_lengths_split() -> None:
    obs = [_obs("AGAT" * 5, name="a"), _obs("AGAT" * 6, name="b")]
    out = cluster_observations(obs)
    assert len(out) == 2
    lengths = sorted(len(c.consensus) for c in out)
    assert lengths == [20, 24]


def test_length_tolerance_merges_close_bins() -> None:
    # Length 48 + length 49 with tolerance 1 → same cluster bin.
    obs = [_obs("AGAT" * 12, name=f"x{i}") for i in range(3)]
    obs.append(_obs("AGAT" * 12 + "A", name="y"))
    out = cluster_observations(obs, len_tolerance_bp=1)
    # Even merged in the same length bin, identity differs (1 bp diff at end),
    # so they split inside the bin unless identity_threshold drops it. With the
    # default 0.97 and 48 vs 49 bp, identity = 1 - 1/49 ≈ 0.98 → merged.
    assert len(out) == 1
    assert out[0].n_reads == 4


def test_zero_length_tolerance_keeps_bins_separate() -> None:
    obs = [
        _obs("AGAT" * 12, name="a"),
        _obs("AGAT" * 12, name="b"),
        _obs("AGAT" * 12 + "A", name="c"),
    ]
    out = cluster_observations(obs, len_tolerance_bp=0)
    lengths = sorted(len(c.consensus) for c in out)
    assert lengths == [48, 49]


def test_identity_threshold_clusters_noise(synth_bam_with_noise: Path) -> None:
    obs = pileup_locus(synth_bam_with_noise, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    out_loose = cluster_observations(obs, identity_threshold=0.95)
    out_strict = cluster_observations(obs, identity_threshold=0.999)
    assert sum(c.n_reads for c in out_loose) == sum(c.n_reads for c in out_strict)
    # With identity 0.95, the single-substitution read joins the main cluster.
    assert len(out_loose) == 1
    # With identity 0.999, the single-substitution read forms its own cluster.
    assert len(out_strict) == 2


def test_clusters_sorted_by_n_reads_desc() -> None:
    obs: list[Observation] = []
    for _ in range(2):
        obs.append(_obs("AGAT" * 5))
    for _ in range(5):
        obs.append(_obs("AGAT" * 6))
    for _ in range(3):
        obs.append(_obs("AGAT" * 7))
    out = cluster_observations(obs)
    assert [c.n_reads for c in out] == [5, 3, 2]


def test_invalid_identity_threshold_rejected() -> None:
    with pytest.raises(EvidenceError):
        cluster_observations([_obs("A")], identity_threshold=0.0)
    with pytest.raises(EvidenceError):
        cluster_observations([_obs("A")], identity_threshold=1.5)


def test_invalid_len_tolerance_rejected() -> None:
    with pytest.raises(EvidenceError):
        cluster_observations([_obs("A")], len_tolerance_bp=-1)


def test_cluster_haplotype_aggregation() -> None:
    obs = [
        _obs("AGAT" * 5, hp=1, name="a1"),
        _obs("AGAT" * 5, hp=1, name="a2"),
        _obs("AGAT" * 5, hp=2, name="a3"),
        _obs("AGAT" * 5, hp=None, name="a4"),
    ]
    out = cluster_observations(obs)
    assert len(out) == 1
    c = out[0]
    assert c.n_reads == 4
    assert c.n_hp1 == 2
    assert c.n_hp2 == 1
    assert c.n_hp_none == 1
    assert c.read_ids == ["a1", "a2", "a3", "a4"]


def test_full_pipeline_pileup_then_cluster(synth_bam_heterozygous: Path) -> None:
    """End-to-end: BAM → pileup → cluster reproduces the synthetic profile."""
    obs = pileup_locus(synth_bam_heterozygous, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    clusters = cluster_observations(obs)
    assert len(clusters) == 3
    # Sorted by n_reads desc: CE12 (5), CE11 (4), CE10 (2)
    assert clusters[0].n_reads == 5
    assert clusters[0].n_hp1 == 5
    assert clusters[0].n_hp2 == 0
    assert clusters[0].consensus == "AGAT" * 12

    assert clusters[1].n_reads == 4
    assert clusters[1].n_hp2 == 4
    assert clusters[1].consensus == "AGAT" * 11

    assert clusters[2].n_reads == 2
    assert clusters[2].n_hp_none == 2
    assert clusters[2].consensus == "AGAT" * 10


def test_cluster_dataclass_is_mutable_for_consensus_polish() -> None:
    """We keep Cluster non-frozen because a future polish step may rewrite consensus."""
    c = Cluster(consensus="A")
    c.consensus = "AGAT"  # must not raise
    assert c.consensus == "AGAT"


def test_cluster_records_how_its_consensus_was_derived() -> None:
    """The record must distinguish a polished consensus from a single read."""
    from frontstr.evidence.cluster import ConsensusMethod

    single = cluster_observations([_obs("AGAT" * 5)])
    assert single[0].consensus_method == ConsensusMethod.SINGLE

    many = cluster_observations([_obs("AGAT" * 5, name=f"r{i}") for i in range(6)])
    assert many[0].consensus_method in (
        ConsensusMethod.POA_ABPOA,
        ConsensusMethod.POA_SPOA,
        ConsensusMethod.MODE,
    )


# ---------------------------------------------------------------------------
# Repeat-core binning: flank errors must not split an allele, microvariants must
# not be merged. See the module docstring in frontstr/evidence/cluster.py.
# ---------------------------------------------------------------------------

_LFLANK = "GCTTCCGAGTGCAGGTCACAGGGAACACAGACTCCATGGTG"
_RFLANK = "AGGGAAATAAGGGAGGAACAGGCCAATGGGAATCACCCCAG"

# TH01: allele 9 is [AATG]9 (36 bp core); allele 9.3 is [AATG]6 ATG [AATG]3
# (39 bp core). Both are 9 repeat units, so a unit-count key would merge them.
_TH01_9 = "AATG" * 9
_TH01_93 = "AATG" * 6 + "ATG" + "AATG" * 3


def _read(core: str, *, lflank: str = _LFLANK, rflank: str = _RFLANK, name: str = "r"):
    return _obs(lflank + core + rflank, name=name)


def test_flank_indel_does_not_split_an_allele() -> None:
    """The dominant ONT error mode: a 1 bp indel in the ±100 bp flank.

    Binned on raw window length these are two alleles; binned on the repeat
    core they are correctly one.
    """
    reads = [
        _read(_TH01_9, name="clean"),
        _read(_TH01_9, lflank=_LFLANK[:-1], name="del_left"),
        _read(_TH01_9, rflank=_RFLANK + "T", name="ins_right"),
    ]
    assert len(cluster_observations(reads)) == 3, "raw-length binning splits them"

    out = cluster_observations(reads, motifs=["AATG"])
    assert len(out) == 1
    assert out[0].n_reads == 3


def test_microvariant_is_not_merged_into_the_full_repeat() -> None:
    """TH01 9 vs 9.3 — same unit count, different core length. Must stay apart."""
    reads = [_read(_TH01_9, name=f"a{i}") for i in range(6)]
    reads += [_read(_TH01_93, name=f"b{i}") for i in range(5)]

    out = cluster_observations(reads, motifs=["AATG"])
    assert len(out) == 2
    assert sorted(c.n_reads for c in out) == [5, 6]


def test_a_real_neighbouring_allele_stays_separate() -> None:
    """9 vs 10 repeat units differ by a whole unit — never merge."""
    reads = [_read("AATG" * 9, name=f"a{i}") for i in range(5)]
    reads += [_read("AATG" * 10, name=f"b{i}") for i in range(5)]
    assert len(cluster_observations(reads, motifs=["AATG"])) == 2


def test_without_motifs_behaviour_is_unchanged() -> None:
    """Callers that pass no motifs keep the original window-length binning."""
    reads = [_read(_TH01_9, name="a"), _read(_TH01_9, lflank=_LFLANK[:-1], name="b")]
    assert len(cluster_observations(reads)) == 2


def test_read_without_a_detectable_core_is_kept() -> None:
    """No motif run must mean 'fall back to window length', never 'drop the read'."""
    reads = [_read(_TH01_9, name="ok"), _obs("GCGCGCGCGCGCGCGCGCGC", name="nocore")]
    out = cluster_observations(reads, motifs=["AATG"])
    assert sum(c.n_reads for c in out) == 2


def test_reverse_strand_core_is_located_on_the_canonical_orientation() -> None:
    """Minus-strand markers must find their core after reverse-complementing."""
    import frontstr.motifs as m

    rc_reads = [
        _obs(m.reverse_complement(_LFLANK + _TH01_9 + _RFLANK), name="a"),
        _obs(m.reverse_complement(_LFLANK[:-1] + _TH01_9 + _RFLANK), name="b"),
    ]
    out = cluster_observations(rc_reads, motifs=["AATG"], strand="-")
    assert len(out) == 1
