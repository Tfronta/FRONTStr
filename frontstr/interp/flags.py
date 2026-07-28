"""Derivation of structured :class:`Flag` records from a finished marker call.

Flags are the machine-readable, auditable replacement for free-text warnings.
Producers append flags at the point a condition is decided:

- :func:`frontstr.interp.concordance.cross_check` emits ``LONGTR_DISCORDANT``.
- :func:`derive_marker_flags` (called at the end of ``interpret_marker``) emits
  the conditions that are intrinsic to the finished :class:`MarkerResult`
  (triallelic / mixture / iso-allele / unpolished consensus). Coverage /
  dropout / strand-bias flags that need a run-level threshold are added by the
  QC layer, not here.

Derivation is idempotent: a code already present on the result is not added
again, so callers may run it after other producers without duplication.
"""

from __future__ import annotations

from frontstr.evidence.consensus import ConsensusMethod
from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    Flag,
    FlagCode,
    MarkerResult,
    TriType,
)


def derive_marker_flags(result: MarkerResult) -> None:
    """Append intrinsic marker-level flags to ``result.flags`` in place."""
    present = {f.code for f in result.flags}

    def add(code: FlagCode, message: str) -> None:
        if code not in present:
            result.flags.append(Flag.of(code, message))
            present.add(code)

    if result.tri_type == TriType.MIXTURE_SUSPECTED:
        add(
            FlagCode.MIXTURE_SUSPECTED,
            "More than two alleles present at a locus without allow_triallelic; "
            "possible mixture — review.",
        )
    elif result.tri_type in (TriType.TYPE_I_UNBALANCED, TriType.TYPE_II_BALANCED):
        add(
            FlagCode.TRIALLELIC,
            f"Triallelic pattern ({result.tri_type.value}) at {result.marker_name}.",
        )

    if _mark_isoalleles(result.alleles_called):
        add(
            FlagCode.ISOALLELE,
            f"Iso-alleles at {result.marker_name}: same allele number, "
            "different sequence — see the Sequences view.",
        )

    phantoms = [a for a in result.alleles if a.status == AlleleStatus.HP_PHANTOM]
    if phantoms:
        detail = ", ".join(f"{a.length_bp} bp / {a.n_reads_total} reads" for a in phantoms)
        add(
            FlagCode.HP_PHANTOM_COLLAPSED,
            f"{len(phantoms)} candidate(s) at {result.marker_name} suppressed as "
            f"same-haplotype splits of a stronger allele ({detail}); their reads "
            "are counted in the owning allele's n_reads_absorbed.",
        )

    split_blocks = [a for a in result.alleles if a.n_phase_sets > 1]
    if split_blocks:
        detail = ", ".join(
            f"cluster {a.cluster_index} ({a.n_phase_sets} blocks)" for a in split_blocks
        )
        add(
            FlagCode.PHASE_BLOCK_SPLIT,
            f"{result.marker_name}: haplotype-tagged reads span more than one phase "
            f"block ({detail}). HP labels are local to a block, so those clusters "
            "were treated as unphased — do not read their HP counts as haplotypes.",
        )

    rescued = [a for a in result.alleles_called if a.hp_rescued]
    if rescued:
        major = result.alleles_called[0]
        detail = ", ".join(
            f"{a.number_label or f'cluster {a.cluster_index}'} on "
            f"{a.n_reads_total} reads vs {major.n_reads_total}"
            for a in rescued
        )
        add(
            FlagCode.HP_RESCUED_HET,
            f"{result.marker_name} is called heterozygous on phasing rather than "
            f"peak balance ({detail}): the allele(s) sit on the opposite haplotype "
            "from the major one, so the read imbalance is coverage, not homozygosity.",
        )

    unpolished = [a for a in result.alleles_called if a.consensus_method == ConsensusMethod.MODE]
    if unpolished:
        add(
            FlagCode.CONSENSUS_FALLBACK,
            f"No POA backend installed: {len(unpolished)} called allele(s) at "
            f"{result.marker_name} use an unpolished single-read consensus. "
            "ISFG string and iso-allele match are not reliable — install "
            "pyabpoa or pyspoa (pip install 'frontstr[poa]').",
        )


def _mark_isoalleles(called: list[Allele]) -> bool:
    """Set ``iso.is_isoallele`` on involved called alleles; return True if any.

    An allele is an iso-allele when it shares its (canonical) number with a
    sibling of different ISFG structure, or the catalog resolved it to a named
    iso-variant (``iso.suffix`` set).
    """
    found = False
    by_number: dict[float, list[Allele]] = {}
    for a in called:
        if a.number is not None:
            by_number.setdefault(a.number, []).append(a)
    for group in by_number.values():
        if len({a.isfg for a in group}) > 1:
            for a in group:
                a.iso.is_isoallele = True
            found = True
    for a in called:
        if a.iso.suffix:
            a.iso.is_isoallele = True
            found = True
    return found
