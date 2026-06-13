"""Derivation of structured :class:`Flag` records from a finished marker call.

Flags are the machine-readable, auditable replacement for free-text warnings.
Producers append flags at the point a condition is decided:

- :func:`frontstr.interp.concordance.cross_check` emits ``LONGTR_DISCORDANT``.
- :func:`derive_marker_flags` (called at the end of ``interpret_marker``) emits
  the conditions that are intrinsic to the finished :class:`MarkerResult`
  (triallelic / mixture). Coverage / dropout / strand-bias flags that need a
  run-level threshold are added by the QC layer, not here.

Derivation is idempotent: a code already present on the result is not added
again, so callers may run it after other producers without duplication.
"""

from __future__ import annotations

from frontstr.interp.models import Allele, Flag, FlagCode, MarkerResult, TriType


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
