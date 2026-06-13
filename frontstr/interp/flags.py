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

from frontstr.interp.models import Flag, FlagCode, MarkerResult, TriType


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
