"""Cluster observations into alleles with ONT-aware tolerance.

Two-stage clustering (plan-longtr-improved.md §5.3):

1. **Length binning** — observations of equal (or near-equal, within
   ``len_tolerance_bp``) length are gathered together.
2. **Identity merge inside the bin** — seed-and-grow over edlib edit distance;
   anything within ``identity_threshold`` of the seed joins the cluster.
3. **Consensus** — POA via ``pyabpoa`` if installed, else the most common
   sequence in the cluster (lossy fallback adequate for HiFi / clean ONT).

Result objects carry per-haplotype tallies so the report can render HP-stacks
without re-iterating the BAM.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field

from frontstr.errors import EvidenceError
from frontstr.evidence.pileup import Observation


@dataclass(slots=True)
class Cluster:
    """A group of observations that represent a single allele candidate."""

    consensus: str
    members: list[Observation] = field(default_factory=list)

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
) -> list[Cluster]:
    """Cluster :class:`Observation` instances by length and sequence identity.

    Args:
        obs: Observations from :func:`pileup_locus`.
        len_tolerance_bp: Merge length bins that differ by at most this many bp.
            Use ``0`` for compound motifs (must not collapse isoalleles); use
            ``1`` or ``2`` for homopolymer-heavy loci on ONT.
        identity_threshold: Pairwise identity (``1 - edit_distance / max_len``)
            required to join an existing cluster's seed.

    Returns:
        Clusters sorted by ``n_reads`` descending. Each consensus string is in
        reference orientation.
    """
    if not 0.0 < identity_threshold <= 1.0:
        raise EvidenceError(
            f"identity_threshold must be in (0, 1], got {identity_threshold!r}"
        )
    if len_tolerance_bp < 0:
        raise EvidenceError(
            f"len_tolerance_bp must be >= 0, got {len_tolerance_bp!r}"
        )

    if not obs:
        return []

    bins: dict[int, list[Observation]] = defaultdict(list)
    for o in obs:
        bins[len(o.sequence)].append(o)

    if len_tolerance_bp > 0:
        bins = _merge_close_length_bins(bins, len_tolerance_bp)

    clusters: list[Cluster] = []
    for members in bins.values():
        clusters.extend(_cluster_by_identity(members, identity_threshold))

    clusters.sort(key=lambda c: c.n_reads, reverse=True)
    return clusters


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
    members: list[Observation], identity_threshold: float
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
        consensus = _consensus_of([m.sequence for m in cluster_members])
        clusters.append(Cluster(consensus=consensus, members=cluster_members))
    return clusters


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


def _consensus_of(seqs: list[str]) -> str:
    """Compute a consensus sequence for a cluster.

    Strategy:

    1. If only one sequence, return it verbatim.
    2. If :mod:`pyabpoa` is installed, run POA and return its consensus.
    3. Otherwise fall back to the most common sequence in the group
       (sufficient for HiFi / clean ONT; only used when POA is unavailable
       e.g. on macOS arm64 without prebuilt wheels).
    """
    if not seqs:
        return ""
    if len(seqs) == 1:
        return seqs[0]

    poa = _poa_aligner()
    if poa is not None:
        try:
            result = poa.msa(seqs, out_cons=True, out_msa=False)
            cons = getattr(result, "cons_seq", None)
            if cons:
                return str(cons[0])
        except Exception:
            pass

    counter = Counter(seqs)
    return counter.most_common(1)[0][0]


_POA_SINGLETON: object | None | bool = False  # False = not yet attempted


def _poa_aligner() -> object | None:
    """Lazy-load the POA aligner, returning ``None`` if pyabpoa is unavailable."""
    global _POA_SINGLETON
    if _POA_SINGLETON is False:
        try:
            import pyabpoa  # type: ignore[import-not-found]

            _POA_SINGLETON = pyabpoa.msa_aligner()
        except ImportError:
            _POA_SINGLETON = None
    return _POA_SINGLETON  # type: ignore[return-value]
