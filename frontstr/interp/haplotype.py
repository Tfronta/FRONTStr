"""Haplotype-aware suppression of split-allele phantoms.

The biological invariant
------------------------

At an autosomal locus a diploid individual carries **exactly one allele per
haplotype**. So when two candidate clusters are both confidently assigned to
the *same* haplotype and differ by less than one repeat unit, at most one of
them can be a real allele — the other is the same allele split in two by
sequencing error.

Two mechanisms produce that split on ONT data, both observed on the R10/Dorado
1000G slices:

1. **Identity split.** Clustering bins by exact length first, then seeds and
   grows by pairwise edit-distance identity. Two *raw reads* of the same allele
   each carrying ~1-2% error differ from each other by ~2-4%, which sits right
   on the 0.97 default threshold — so a noisy seed read rejects its own
   siblings. Signature: identical length, identical dominant haplotype
   (HG00097 D8S1179, HG00263 D8S1179/D16S539, HG00097 D13S317).
2. **Length split.** A single-bp indel error moves a read into the adjacent
   length bin, which ``len_tolerance_bp=0`` never re-merges. Signature: 1 bp
   apart, identical dominant haplotype (GM19038 D12S391 289 bp / 288 bp, both
   HP1).

Why this is the right fix rather than a read-count floor
--------------------------------------------------------

A flat "a 3rd allele needs ≥ N reads" floor cannot separate the two cases: in
GM19038 D12S391 the two strongest clusters are *both* HP1 phantoms of one
allele and the genuine HP2 allele is the 5-read cluster. Raising the floor
would delete the real allele and emit a confidently wrong homozygote. The
haplotype tag is orthogonal evidence that resolves it correctly.

Conservatism
------------

This suppression only ever fires when all of the following hold, so it cannot
mask a genuine second contributor or a genuine triallelic pattern:

- both clusters carry at least ``min_tagged_reads`` phased reads and have a
  *confident* dominant haplotype (``hp_purity``) — on an unphased BAM no
  cluster qualifies and the whole step is a no-op;
- they share that haplotype;
- they differ by at most ``max_phantom_bp`` — well under one repeat unit. A
  second contributor or a true third allele differs by a whole repeat unit.
- the marker does not set ``allow_triallelic`` (a duplication genuinely puts
  two alleles on one haplotype).

Suppressed clusters are demoted to :class:`AlleleStatus.HP_PHANTOM` rather
than deleted: they stay in the record, the owning allele records how many
reads were absorbed, and the marker carries a ``HP_PHANTOM_COLLAPSED`` flag so
the decision is auditable.
"""

from __future__ import annotations

from frontstr.interp.models import Allele, AlleleStatus
from frontstr.panel.models import System

#: Max length difference for two same-haplotype clusters to be the same allele.
#: Kept far below one repeat unit (4 bp for a tetramer) so a genuine adjacent
#: allele is never absorbed. Observed phantoms are 0-1 bp apart.
DEFAULT_MAX_PHANTOM_BP = 2

#: Fraction of a cluster's *tagged* reads that must sit on one haplotype before
#: we treat that cluster as belonging to it.
DEFAULT_HP_PURITY = 0.8

#: Tagged reads a cluster needs before its haplotype assignment is evidence.
#: This is deliberately a per-cluster count, not a locus-wide tagged fraction:
#: whatshap leaves a read untagged when it cannot place it, which is missing
#: evidence, not evidence against. Gating on the locus fraction rejected
#: unambiguous cases (HG00097 D13S317: two 259 bp clusters, 100% HP1 each, but
#: only 44% of the locus tagged).
DEFAULT_MIN_TAGGED_READS = 3


def dominant_hp(
    allele: Allele,
    *,
    hp_purity: float = DEFAULT_HP_PURITY,
    min_tagged_reads: int = DEFAULT_MIN_TAGGED_READS,
) -> int | None:
    """Return ``1`` / ``2`` when the cluster's tagged reads agree, else ``None``.

    Untagged reads are ignored rather than counted against purity: a cluster of
    9 HP1 reads plus 3 untagged reads is still an HP1 cluster. Clusters with
    fewer than ``min_tagged_reads`` tagged reads return ``None`` — too little
    evidence to claim a haplotype, which is also what makes this a no-op on
    unphased BAMs.
    """
    tagged = allele.n_reads_hp1 + allele.n_reads_hp2
    if tagged < min_tagged_reads:
        return None
    if allele.n_reads_hp1 / tagged >= hp_purity:
        return 1
    if allele.n_reads_hp2 / tagged >= hp_purity:
        return 2
    return None


def suppress_hp_phantoms(
    alleles: list[Allele],
    system: System,
    *,
    max_phantom_bp: int = DEFAULT_MAX_PHANTOM_BP,
    hp_purity: float = DEFAULT_HP_PURITY,
    min_tagged_reads: int = DEFAULT_MIN_TAGGED_READS,
) -> int:
    """Demote same-haplotype split-allele phantoms in place.

    For each haplotype, the strongest candidate is the owner; any weaker
    candidate on the same haplotype within ``max_phantom_bp`` of the owner is
    demoted to :class:`AlleleStatus.HP_PHANTOM` and its reads are recorded on
    the owner's ``n_reads_absorbed``.

    Args:
        alleles: All classified candidates for the marker; mutated in place.
        system: Marker definition. Markers with ``allow_triallelic`` are
            skipped — a duplication legitimately puts two alleles on one
            haplotype.
        max_phantom_bp: Max length difference between owner and phantom.
        hp_purity: Fraction of tagged reads needed to assign a cluster to a
            haplotype.
        min_tagged_reads: Tagged reads a cluster needs before its haplotype is
            treated as evidence.

    Returns:
        Number of candidates demoted. ``0`` means nothing changed, which is the
        case for unphased BAMs and for triallelic-capable markers.
    """
    if system.allow_triallelic:
        return 0

    candidates = [
        a
        for a in alleles
        if a.status in (AlleleStatus.ALLELE, AlleleStatus.INEXACT_ALLELE) and not a.is_deletion
    ]
    if len(candidates) < 2:
        return 0

    by_hp: dict[int, list[Allele]] = {}
    for a in candidates:
        hp = dominant_hp(a, hp_purity=hp_purity, min_tagged_reads=min_tagged_reads)
        if hp is not None:
            by_hp.setdefault(hp, []).append(a)

    demoted = 0
    for group in by_hp.values():
        group.sort(key=lambda a: (-a.n_reads_total, a.cluster_index))
        owner = group[0]
        for phantom in group[1:]:
            if abs(phantom.length_bp - owner.length_bp) > max_phantom_bp:
                continue  # a whole repeat unit apart — a different allele
            phantom.status = AlleleStatus.HP_PHANTOM
            owner.n_reads_absorbed += phantom.n_reads_total
            demoted += 1
    return demoted
