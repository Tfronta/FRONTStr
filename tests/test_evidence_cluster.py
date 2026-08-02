"""Tests for :func:`frontstr.evidence.cluster.cluster_observations`."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.errors import EvidenceError
from frontstr.evidence.cluster import (
    Cluster,
    _merge_identical_consensus,
    cluster_observations,
)
from frontstr.evidence.consensus import poa_backend_name
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


def _obs_at(seq: str, n: int, *, start: int = 0, hp: int | None = None) -> list[Observation]:
    return [
        Observation(
            read_id=f"r{start + i}",
            sequence=seq,
            hp=hp,
            ps=1 if hp else None,
            mean_qual=30.0,
            strand="+",
            flank_left_ok=True,
            flank_right_ok=True,
        )
        for i in range(n)
    ]


def test_identical_consensus_clusters_are_folded_together() -> None:
    """One allele that seeded two clusters is one allele again after POA.

    Seed-and-grow compares raw read to raw read, and two ONT reads of the same
    allele differ by 2-4%, so one allele can start two clusters. POA polishes
    each separately and both land on the same sequence. Byte-identical
    consensus is not a threshold question: it is the same allele.
    """
    seq = "A" * 60 + "AATG" * 12 + "T" * 60
    clusters = cluster_observations(_obs_at(seq, 8), motifs=["AATG"])
    assert len(clusters) == 1
    assert clusters[0].n_reads == 8

    # Force the split the identity stage would produce, then check the merge
    # puts the reads back on one cluster rather than leaving two.
    merged: dict[str, int] = {}
    split = [
        Cluster(consensus=seq, members=_obs_at(seq, 3)),
        Cluster(consensus=seq, members=_obs_at(seq, 3, start=100)),
        Cluster(consensus="A" * 60 + "AATG" * 10 + "T" * 60, members=_obs_at("x", 5)),
    ]
    out = _merge_identical_consensus(split, merged)
    assert len(out) == 2
    assert out[0].n_reads == 6, "the two fragments must add up"
    assert merged == {seq: 2}


def test_merge_keeps_genuinely_different_alleles_apart() -> None:
    """A real heterozygote must survive the merge untouched."""
    a = "A" * 60 + "AATG" * 12 + "T" * 60
    b = "A" * 60 + "AATG" * 9 + "T" * 60
    out = _merge_identical_consensus(
        [Cluster(consensus=a, members=_obs_at(a, 10)), Cluster(consensus=b, members=_obs_at(b, 9))]
    )
    assert len(out) == 2
    assert {c.consensus for c in out} == {a, b}


def test_merge_leaves_empty_consensus_clusters_alone() -> None:
    """Two empty consensuses are two failures, not one shared allele."""
    out = _merge_identical_consensus(
        [
            Cluster(consensus="", members=_obs_at("x", 1)),
            Cluster(consensus="", members=_obs_at("y", 1)),
        ]
    )
    assert len(out) == 2


def test_merged_fragments_carry_their_haplotype_reads_across() -> None:
    """The point of the merge on the cohort was the haplotype tallies.

    HG04161 FGA split one allele into two 3-read clusters. Apart, neither
    cleared the heterozygote ratio against a 13-read major and the locus was
    called homozygous; together they are 6 of 13, a ratio of 0.46. The merged
    cluster has to carry both fragments' HP counts or the rescue still cannot
    see them.
    """
    seq = "A" * 60 + "AAAG" * 24 + "T" * 60
    out = _merge_identical_consensus(
        [
            Cluster(consensus=seq, members=_obs_at(seq, 3, hp=2)),
            Cluster(consensus=seq, members=_obs_at(seq, 3, start=50, hp=2)),
        ]
    )
    assert len(out) == 1
    assert out[0].n_reads == 6
    assert out[0].n_hp2 == 6


# --- Refinement against the consensus -------------------------------------

poa_only = pytest.mark.skipif(
    not poa_backend_name(), reason="no POA backend installed (pyabpoa / pyspoa)"
)


def _noisy(truth: str, positions: list[int]) -> str:
    """``truth`` with a substitution at each position, length preserved.

    Substitutions only, so every read stays in one length bin and the test is
    about identity rather than binning.
    """
    out = list(truth)
    for p in positions:
        out[p] = "C" if out[p] != "C" else "G"
    return "".join(out)


_TRUTH = "AATG" * 50  # 200 bp; at 0.97 identity the budget is 6 edits.


def test_a_read_moves_to_the_consensus_it_actually_matches() -> None:
    """The decision rule, isolated from how good POA's consensus is."""
    from frontstr.evidence.cluster import _refine_by_consensus

    stray = _obs(_noisy(_TRUTH, [7]), name="stray")
    home = _obs_at(_TRUTH, 3)
    clusters = [
        Cluster(consensus=stray.sequence, members=[stray]),
        Cluster(consensus=_TRUTH, members=home),
    ]
    out, moved = _refine_by_consensus(clusters, [stray, *home], 0.97)

    assert moved == 1
    assert len(out) == 1, "the stray read left no cluster behind"
    assert out[0].n_reads == 4
    assert "stray" in out[0].read_ids


def test_a_genuinely_different_sequence_keeps_its_cluster() -> None:
    """Refinement must not become a merge-everything pass."""
    from frontstr.evidence.cluster import _refine_by_consensus

    other = _obs(_noisy(_TRUTH, list(range(0, 60, 4))), name="other")
    home = _obs_at(_TRUTH, 3)
    clusters = [
        Cluster(consensus=other.sequence, members=[other]),
        Cluster(consensus=_TRUTH, members=home),
    ]
    out, moved = _refine_by_consensus(clusters, [other, *home], 0.97)

    assert moved == 0
    assert len(out) == 2


def test_a_single_read_never_attracts() -> None:
    """A lone read's "consensus" is itself, so it is not evidence of anything.

    Were singletons allowed to attract, every read would match its own cluster
    perfectly and nothing could ever move.
    """
    from frontstr.evidence.cluster import _refine_by_consensus

    a, b = _obs(_noisy(_TRUTH, [3]), name="a"), _obs(_noisy(_TRUTH, [9]), name="b")
    out, moved = _refine_by_consensus(
        [Cluster(consensus=a.sequence, members=[a]), Cluster(consensus=b.sequence, members=[b])],
        [a, b],
        0.97,
    )

    assert moved == 0
    assert len(out) == 2


@poa_only
def test_a_noisy_seed_no_longer_splits_one_allele() -> None:
    """The defect: which reads cluster depended on which read seeded.

    Every read here is the same allele. The first one carries five errors, so
    seed-and-grow measures the rest against those errors as well as their own
    and pushes them past the threshold. The POA consensus of the others does
    not carry them, and against it the seed is back inside.
    """
    reads = [
        _obs(_noisy(_TRUTH, [10, 30, 50, 70, 90]), name="seed"),
        _obs(_noisy(_TRUTH, [11, 31]), name="r1"),
        _obs(_noisy(_TRUTH, [13, 33]), name="r2"),
        _obs(_noisy(_TRUTH, [15, 35]), name="r3"),
        _obs(_noisy(_TRUTH, [17, 37]), name="r4"),
    ]
    out = cluster_observations(reads)

    assert len(out) == 1, "one allele, one cluster"
    assert out[0].n_reads == 5
    assert out[0].consensus == _TRUTH


@poa_only
def test_clustering_no_longer_depends_on_read_order() -> None:
    """Checkable without a comparator, which is the point.

    The same reads arriving in a different order are the same sample, so they
    have to yield the same clusters.
    """
    reads = [
        _obs(_noisy(_TRUTH, [10, 30, 50, 70, 90]), name="seed"),
        _obs(_noisy(_TRUTH, [11, 31]), name="r1"),
        _obs(_noisy(_TRUTH, [13, 33]), name="r2"),
        _obs(_noisy(_TRUTH, [15, 35]), name="r3"),
        _obs(_noisy(_TRUTH, [17, 37]), name="r4"),
    ]
    seed_first = cluster_observations(reads)
    seed_last = cluster_observations([*reads[1:], reads[0]])

    assert [c.n_reads for c in seed_first] == [c.n_reads for c in seed_last]
    assert [c.consensus for c in seed_first] == [c.consensus for c in seed_last]
    assert sorted(seed_first[0].read_ids) == sorted(seed_last[0].read_ids)


@poa_only
def test_two_alleles_stay_two_alleles_under_refinement() -> None:
    """Refinement runs inside a length bin, so a real heterozygote is safe.

    Noise goes in the flank, where ONT actually puts most of it and where it
    leaves the repeat core, and therefore the binning key, untouched.
    """
    flank_l, flank_r = "GCTA" * 25, "TCGA" * 25
    short = flank_l + "AATG" * 10 + flank_r
    long_ = flank_l + "AATG" * 14 + flank_r
    reads = [
        *[_obs(_noisy(short, [2 * i + 1]), name=f"s{i}") for i in range(4)],
        *[_obs(_noisy(long_, [2 * i + 2]), name=f"l{i}") for i in range(4)],
    ]
    out = cluster_observations(reads, motifs=["AATG"])

    assert len(out) == 2
    assert sorted(c.n_reads for c in out) == [4, 4]
    assert {c.consensus for c in out} == {short, long_}


@poa_only
def test_refinement_reports_every_read_it_moved() -> None:
    """The trace has to be able to say a read changed candidate allele."""
    reads = [
        _obs(_noisy(_TRUTH, [10, 30, 50, 70, 90]), name="seed"),
        *[_obs(_noisy(_TRUTH, [11 + 2 * i, 31 + 2 * i]), name=f"r{i}") for i in range(4)],
    ]
    reassigned: dict[int, int] = {}
    cluster_observations(reads, reassigned=reassigned)

    assert sum(reassigned.values()) == 1
    assert list(reassigned) == [len(_TRUTH)], "keyed by the bin the move happened in"
