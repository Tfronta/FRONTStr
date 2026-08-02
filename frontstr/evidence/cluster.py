"""Cluster observations into alleles with ONT-aware tolerance.

Clustering runs in five stages:

1. **Length binning** — observations of equal (or near-equal, within
   ``len_tolerance_bp``) length are gathered together. When ``motifs`` are
   supplied the binning key is the **repeat-core length** rather than the raw
   window length; see :func:`_binning_key`.
2. **Identity merge inside the bin** — seed-and-grow over edlib edit distance;
   anything within ``identity_threshold`` of the seed joins the cluster.
3. **Consensus** — partial-order alignment via
   :mod:`frontstr.evidence.consensus`, which records *how* the consensus was
   derived on each cluster (see :class:`ConsensusMethod`).
4. **Refine against the consensus** — stage 2 measured every read against a raw
   seed read, which carries that read's own errors. Once a consensus exists it
   is the better yardstick, so reads are reassigned against it and the
   consensus rebuilt. See :func:`_refine_by_consensus`.
5. **Merge on identical consensus** — one allele can still seed two clusters;
   once POA has polished both, an identical consensus is proof they were one
   allele all along. See :func:`_merge_identical_consensus`.

Result objects carry per-haplotype tallies so the report can render HP-stacks
without re-iterating the BAM.

Why the binning key is the repeat core
--------------------------------------

Panel windows are the repeat array plus ±100 bp of flank, so roughly 80% of a
window is flank. Binning on raw window length makes every ONT indel error in
that flank — the large majority of errors — split one real allele into two
clusters, which then surfaces downstream as a phantom third allele. Binning on
the repeat-core length ignores flank errors entirely while staying exact inside
the array, so genuine microvariants survive: ``[AATG]9`` (36 bp core) and
``[AATG]6 ATG [AATG]3`` (39 bp core) still land in different bins even though
both are 9 repeat units. That distinction is TH01 9 vs 9.3 — precisely the
thing a length-rounding scheme would destroy.
"""

from __future__ import annotations

import time
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field

from frontstr.errors import EvidenceError
from frontstr.evidence.consensus import ConsensusMethod, build_consensus
from frontstr.evidence.pileup import Observation
from frontstr.motifs import repeat_core_length


@dataclass(slots=True)
class Cluster:
    """A group of observations that represent a single allele candidate."""

    consensus: str
    members: list[Observation] = field(default_factory=list)
    #: How :attr:`consensus` was derived. ``MODE`` means no POA backend was
    #: installed and the consensus is an unpolished single read — the
    #: interpretation layer raises ``CONSENSUS_FALLBACK`` for those.
    consensus_method: ConsensusMethod = ConsensusMethod.SINGLE

    @property
    def n_reads(self) -> int:
        return len(self.members)

    @property
    def n_hp1(self) -> int:
        return sum(1 for o in self.members if o.hp == 1)

    @property
    def n_hp2(self) -> int:
        return sum(1 for o in self.members if o.hp == 2)

    @property
    def n_hp_none(self) -> int:
        return sum(1 for o in self.members if o.hp is None)

    @property
    def phase_sets(self) -> Counter[int]:
        """Phase blocks the haplotype-tagged members came from, with counts.

        More than one entry means this cluster's ``HP`` labels are not
        comparable with each other: HP2 in one block and HP2 in the next are
        unrelated labels. The interpretation layer treats such a cluster as
        unphased rather than trusting the majority.
        """
        return Counter(o.ps for o in self.members if o.hp is not None and o.ps is not None)

    @property
    def n_forward(self) -> int:
        return sum(1 for o in self.members if o.strand == "+")

    @property
    def n_reverse(self) -> int:
        return sum(1 for o in self.members if o.strand == "-")

    @property
    def mean_qual(self) -> float:
        if not self.members:
            return 0.0
        return sum(o.mean_qual for o in self.members) / len(self.members)

    @property
    def read_ids(self) -> list[str]:
        return [o.read_id for o in self.members]


_DEFAULT_IDENTITY_THRESHOLD = 0.97  # ONT R10 simplex; use 0.99 for HiFi / duplex.


def cluster_observations(
    obs: list[Observation],
    *,
    len_tolerance_bp: int = 0,
    identity_threshold: float = _DEFAULT_IDENTITY_THRESHOLD,
    motifs: Sequence[str] | None = None,
    strand: str = "+",
    bin_sizes: dict[int, int] | None = None,
    core_binned: dict[int, int] | None = None,
    consensus_seconds: list[float] | None = None,
    merged_consensus: dict[str, int] | None = None,
    reassigned: dict[int, int] | None = None,
    split_sizes: dict[int, int] | None = None,
) -> list[Cluster]:
    """Cluster :class:`Observation` instances by length and sequence identity.

    Args:
        obs: Observations from :func:`pileup_locus`.
        len_tolerance_bp: Merge length bins that differ by at most this many bp.
            Leave at ``0``: with repeat-core binning the flank errors this was
            meant to absorb are already gone, and a non-zero tolerance would
            merge genuine microvariants (TH01 9 vs 9.3 are 3 bp apart).
        identity_threshold: Pairwise identity (``1 - edit_distance / max_len``)
            required to join an existing cluster's seed.
        motifs: Marker motifs. When given, reads are binned by repeat-core
            length instead of raw window length, so indel errors in the ±100 bp
            flanks stop splitting real alleles. Reads with no detectable core
            fall back to their window length.
        strand: ``"-"`` for markers whose canonical motif reads on the minus
            strand; the sequence is reverse-complemented before the core is
            located.
        bin_sizes: Optional dict filled with ``{binning key: n reads}``, so a
            trace can show why reads grouped as they did.
        consensus_seconds: Optional list appended with the wall-clock seconds
            each POA consensus took, so a trace can attribute clustering time
            between binning and consensus. ``None`` skips the clock entirely.
        core_binned: Optional dict filled with ``{binning key: n reads whose
            repeat core was actually locatable}``. A key whose count is lower
            than ``bin_sizes`` contains reads binned on raw window length, which
            includes flank error — worth surfacing rather than hiding.
        merged_consensus: Optional dict filled with ``{consensus: n clusters}``
            for each set of fragments folded together by
            :func:`_merge_identical_consensus`, so a trace can say a merge
            happened instead of just reporting fewer clusters than bins.
        reassigned: Optional dict filled with ``{binning key: n reads that
            changed cluster}`` during :func:`_refine_by_consensus`. A read
            moving between candidate alleles is exactly the kind of step the
            audit trail has to show rather than absorb.
        split_sizes: Optional dict filled with ``{binning key: n clusters}`` as
            the identity split left them, *before* refinement and the merge.
            The returned list has been through both, so a trace that wants to
            report what splitting alone did cannot count it from the result.

    Returns:
        Clusters sorted by ``n_reads`` descending. Each consensus string is in
        reference orientation.
    """
    if not 0.0 < identity_threshold <= 1.0:
        raise EvidenceError(f"identity_threshold must be in (0, 1], got {identity_threshold!r}")
    if len_tolerance_bp < 0:
        raise EvidenceError(f"len_tolerance_bp must be >= 0, got {len_tolerance_bp!r}")

    if not obs:
        return []

    motif_list = [m for m in (motifs or []) if m]
    bins: dict[int, list[Observation]] = defaultdict(list)
    for o in obs:
        key, from_core = _binning_key(o.sequence, motif_list, strand)
        bins[key].append(o)
        if core_binned is not None and from_core:
            core_binned[key] = core_binned.get(key, 0) + 1
    if bin_sizes is not None:
        bin_sizes.update({k: len(v) for k, v in bins.items()})

    if len_tolerance_bp > 0:
        bins = _merge_close_length_bins(bins, len_tolerance_bp)

    clusters: list[Cluster] = []
    for key, members in bins.items():
        binned = _cluster_by_identity(members, identity_threshold, consensus_seconds)
        if split_sizes is not None:
            split_sizes[key] = len(binned)
        binned, moved = _refine_by_consensus(binned, members, identity_threshold, consensus_seconds)
        if reassigned is not None and moved:
            reassigned[key] = moved
        clusters.extend(binned)

    clusters = _merge_identical_consensus(clusters, merged_consensus)
    clusters.sort(key=lambda c: c.n_reads, reverse=True)
    return clusters


#: Reassignment passes before the result is accepted as-is. Each pass moves
#: reads and rebuilds the consensus they were compared against, so a later pass
#: can move more. The loop stops as soon as a pass moves nothing, which is the
#: normal exit; this cap only guarantees termination, since membership can in
#: principle alternate between two equally good consensuses. It is not a tuned
#: quantity and no call should depend on its value.
_MAX_REFINE_PASSES = 3


def _refine_by_consensus(
    clusters: list[Cluster],
    bin_members: list[Observation],
    identity_threshold: float,
    consensus_seconds: list[float] | None = None,
) -> tuple[list[Cluster], int]:
    """Reassign reads against polished consensuses instead of raw seed reads.

    Seed-and-grow asks "is this read within ``identity_threshold`` of *that
    read*", and the seed is simply whichever read the BAM happened to yield
    first. An ONT read differs from the truth by 2-4%, so the seed spends part
    of the budget on its own errors and reads of the same allele are pushed
    outside it. Which reads those are depends on which read seeded — the split
    is an artifact of input order, not of the sample.

    A POA consensus over several reads has had those random errors averaged
    out, so read-to-consensus is a strictly better estimate of "same allele?"
    than read-to-read. This reassigns every read to the consensus it matches
    best and rebuilds the consensus from the new membership, which is checkable
    without any comparator: the same reads in a different order now converge on
    the same clusters.

    Only clusters of **two or more reads** attract. A singleton's "consensus"
    is the read itself, so it would match itself perfectly and never move, and
    a lone read has none of the error-averaging that justifies the comparison
    in the first place.

    A read that matches no consensus at ``identity_threshold`` stays where it
    is. Nothing is orphaned and no cluster is created here; genuinely distinct
    sequence keeps its own cluster, and two fragments that converge on the same
    sequence are left for :func:`_merge_identical_consensus` to fold together,
    which needs no threshold at all.

    Args:
        clusters: Clusters from one bin.
        bin_members: Every observation in that bin, in BAM order, so a cluster
            that gains reads can be put back into that order before its
            consensus is rebuilt. POA output depends on input order, and the
            reproducibility guarantee rests on it.
        identity_threshold: Same threshold stage 2 used, now applied to a
            cleaner reference sequence. Reused rather than given its own knob:
            it is the same question asked of better evidence, and against a
            consensus the test is the more conservative of the two.
        consensus_seconds: Appended with the seconds each rebuild took, so the
            trace attributes refinement time rather than hiding it in binning.

    Returns:
        The refined clusters, and how many reads changed cluster in total.
        Members stay in BAM order so the consensus stays reproducible.
    """
    if len(clusters) < 2:
        return clusters, 0

    bam_order = {id(o): i for i, o in enumerate(bin_members)}
    moved_total = 0
    for _ in range(_MAX_REFINE_PASSES):
        attractors = sorted(
            (c for c in clusters if c.n_reads >= 2 and c.consensus),
            key=lambda c: c.n_reads,
            reverse=True,
        )
        if not attractors:
            break

        assigned: dict[int, list[Observation]] = {id(c): [] for c in clusters}
        moved = 0
        for cluster in clusters:
            for o in cluster.members:
                # Ties go to the larger cluster, whose consensus rests on more
                # reads: `attractors` is sorted and only a strict improvement
                # displaces the incumbent, so the outcome is deterministic.
                best, best_identity = attractors[0], _identity(attractors[0].consensus, o.sequence)
                for a in attractors[1:]:
                    identity = _identity(a.consensus, o.sequence)
                    if identity > best_identity:
                        best, best_identity = a, identity
                target = best if best_identity >= identity_threshold else cluster
                if target is not cluster:
                    moved += 1
                assigned[id(target)].append(o)
        if not moved:
            break
        moved_total += moved

        rebuilt: list[Cluster] = []
        for cluster in clusters:
            members = assigned[id(cluster)]
            if not members:
                continue
            if len(members) == cluster.n_reads and all(
                a is b for a, b in zip(members, cluster.members, strict=True)
            ):
                rebuilt.append(cluster)
                continue
            members.sort(key=lambda o: bam_order[id(o)])
            consensus, method = _build_consensus_timed(members, consensus_seconds)
            rebuilt.append(Cluster(consensus=consensus, members=members, consensus_method=method))
        clusters = rebuilt

    return clusters, moved_total


def _merge_identical_consensus(
    clusters: list[Cluster], merged: dict[str, int] | None = None
) -> list[Cluster]:
    """Fold together clusters whose polished consensus is byte-identical.

    Seed-and-grow compares a raw read to a raw read. Two ONT reads of the same
    allele differ from each other by 2-4%, straddling ``identity_threshold``, so
    one allele can start two clusters. POA then polishes each of them
    separately and both converge on the same true sequence. What comes out is
    two clusters with the same consensus, which is not a threshold question:
    two clusters carrying byte-identical sequence are the same allele, by
    definition of allele.

    Measured across the 108-sample 1000G ONT cohort, three of the loci
    FRONTStr got wrong were exactly this, and the symptom went both ways. When
    both fragments were strong they were both called and the genotype carried
    the same number twice: HG02187 TPOX ``10/12/12`` and HG03667 FGA
    ``21/22/22``. When both were weak neither cleared the heterozygote ratio
    and the locus was called homozygous: HG04161 FGA split one allele into two
    3-read clusters that together are 6 of 13 reads, a ratio of 0.46 that
    passes comfortably.

    This does not close the ``identity_threshold`` question. It removes the
    case that needs no threshold to decide; near-identical fragments still
    want consensus-refined reassignment.

    Args:
        clusters: Clusters straight out of :func:`_cluster_by_identity`.
        merged: Optional dict filled with ``{consensus: n clusters folded in}``
            for the fragments that were absorbed, so a trace can report a merge
            rather than silently showing fewer clusters than bins.

    Returns:
        One cluster per distinct consensus, members concatenated, in first-seen
        order. Empty consensuses are left alone: they are deletions or failures
        and carry no evidence that they are the same event.
    """
    by_consensus: dict[str, Cluster] = {}
    out: list[Cluster] = []
    for cluster in clusters:
        if not cluster.consensus:
            out.append(cluster)
            continue
        first = by_consensus.get(cluster.consensus)
        if first is None:
            by_consensus[cluster.consensus] = cluster
            out.append(cluster)
            continue
        first.members.extend(cluster.members)
        if merged is not None:
            merged[cluster.consensus] = merged.get(cluster.consensus, 1) + 1
    return out


def _binning_key(sequence: str, motifs: list[str], strand: str) -> tuple[int, bool]:
    """Length used to bin a read, and whether it came from the repeat core.

    A read whose core cannot be located (heavily corrupted array, or no motifs
    configured) falls back to raw window length rather than being dropped or
    lumped together with unrelated reads. The flag is returned so a trace can
    distinguish the two — a window-length key carries flank error, a core key
    does not.
    """
    if not motifs:
        return len(sequence), False
    core = repeat_core_length(sequence, motifs, strand=strand)
    return (core, True) if core is not None else (len(sequence), False)


def _merge_close_length_bins(
    bins: dict[int, list[Observation]], len_tolerance_bp: int
) -> dict[int, list[Observation]]:
    """Greedy left-to-right merge of length bins within ``len_tolerance_bp``."""
    merged: dict[int, list[Observation]] = {}
    keys = sorted(bins.keys())
    i = 0
    while i < len(keys):
        anchor = keys[i]
        group = list(bins[anchor])
        j = i + 1
        while j < len(keys) and keys[j] - anchor <= len_tolerance_bp:
            group.extend(bins[keys[j]])
            j += 1
        merged[anchor] = group
        i = j
    return merged


def _cluster_by_identity(
    members: list[Observation],
    identity_threshold: float,
    consensus_seconds: list[float] | None = None,
) -> list[Cluster]:
    """Seed-and-grow clustering: each seed is the first uncovered observation."""
    clusters: list[Cluster] = []
    remaining = list(members)
    while remaining:
        seed = remaining.pop(0)
        cluster_members: list[Observation] = [seed]
        leftover: list[Observation] = []
        for m in remaining:
            if _identity(seed.sequence, m.sequence) >= identity_threshold:
                cluster_members.append(m)
            else:
                leftover.append(m)
        remaining = leftover
        consensus, method = _build_consensus_timed(cluster_members, consensus_seconds)
        clusters.append(
            Cluster(
                consensus=consensus,
                members=cluster_members,
                consensus_method=method,
            )
        )
    return clusters


def _build_consensus_timed(
    members: list[Observation], consensus_seconds: list[float] | None
) -> tuple[str, ConsensusMethod]:
    """Build a cluster consensus, charging the time to the consensus stage.

    Timed separately so a trace can tell "many divergent candidates to align"
    apart from "many reads to bin" — they look the same from outside clustering
    and mean very different things.
    """
    if consensus_seconds is None:
        return build_consensus([m.sequence for m in members])
    started = time.perf_counter()
    result = build_consensus([m.sequence for m in members])
    consensus_seconds.append(time.perf_counter() - started)
    return result


def _identity(a: str, b: str) -> float:
    """Levenshtein-based identity, ``1 - edit / max(len(a), len(b))``."""
    if not a or not b:
        return 0.0
    if a == b:
        return 1.0
    try:
        import edlib
    except ImportError as exc:
        raise EvidenceError("edlib is required for cluster identity") from exc
    res = edlib.align(a, b, task="distance", mode="NW")
    d = int(res["editDistance"])
    if d < 0:
        return 0.0
    return 1.0 - d / max(len(a), len(b))


__all__ = ["Cluster", "ConsensusMethod", "cluster_observations"]
