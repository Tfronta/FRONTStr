"""Stutter expectation model (LUS / SLUS) ported from toaSTR v1.

For each parent allele observed in the reads we build a set of **virtual
stutter sequences** by manipulating the longest and second-longest
uninterrupted stretches of the motif (LUS and SLUS). Each virtual stutter
contributes an *expected coverage* that scales as ``ST^k * parent_coverage``
where ``k`` is the number of stutter steps and ``ST`` is the per-step
stutter rate.

We emit -1, -2 and +1 stutters; isometric variants (substitutions inside
the LUS) are skipped because LongTR's INEXACT_ALLELE flag is a stronger
forensic signal.

The function returns a mapping ``virtual_sequence -> expected_coverage``.
When a real cluster's consensus matches one of these keys, its expected
stutter is the value at that key (see :mod:`frontstr.interp.classify`).

Reference: research-toastr.md §8 step 4 and plan-longtr-improved.md §6.2.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from frontstr.evidence.cluster import Cluster
from frontstr.motifs import MotifRun, find_motif_runs, top_two_runs
from frontstr.panel.models import System

__all__ = [
    "MotifRun",
    "build_expected_stutter",
    "find_motif_runs",
    "top_two_runs",
    "virtual_stutters",
]


def _modify_run(sequence: str, run: MotifRun, delta_copies: int) -> str | None:
    """Return ``sequence`` with ``run`` shortened/lengthened by ``delta_copies``.

    Returns ``None`` if the modification would make the run vanish or go
    negative (i.e. ``run.n_copies + delta_copies < 1``).
    """
    new_copies = run.n_copies + delta_copies
    if new_copies < 1:
        return None
    prefix = sequence[: run.start]
    suffix = sequence[run.end :]
    return prefix + run.motif * new_copies + suffix


def virtual_stutters(parent_seq: str, motifs: Iterable[str]) -> Iterator[tuple[str, int, str]]:
    """Yield ``(variant_seq, step, source)`` for each virtual stutter from ``parent_seq``.

    ``step`` is 1 or 2 (forensic stutter levels we model). ``source`` is
    ``"LUS"`` or ``"SLUS"``. Plus-stutters (+1) come back as ``step=1``
    too; the difference between minus and plus contribution lives in the
    coefficient table in :func:`build_expected_stutter`.
    """
    runs = find_motif_runs(parent_seq, motifs)
    lus, slus = top_two_runs(runs)
    for run, source in ((lus, "LUS"), (slus, "SLUS")):
        if run is None or run.n_copies < 2:
            continue
        for delta in (-1, -2, +1):
            variant = _modify_run(parent_seq, run, delta)
            if variant is None or variant == parent_seq:
                continue
            yield variant, abs(delta), source


def build_expected_stutter(
    parents: list[Cluster],
    system: System,
    *,
    default_lus: float = 0.10,
    default_slus: float = 0.05,
    plus_factor: float = 0.5,
) -> dict[str, float]:
    """Compute expected stutter coverage per virtual stutter sequence.

    Args:
        parents: Strong cluster candidates (above ``calling_thresh``).
        system: The marker definition. ``system.stutter_overrides`` may
            contain ``"lus"`` / ``"slus"`` / ``"plus_factor"`` to tune the
            model per-locus.
        default_lus: Per-step rate inside the LUS (default 10%).
        default_slus: Per-step rate inside the SLUS (default 5%).
        plus_factor: Scale for ``+1`` stutters relative to ``-1`` (default 0.5,
            reflecting that forward stutters are roughly half as frequent).

    Returns:
        Mapping ``virtual_sequence -> expected_coverage``. The same key may
        be hit by multiple parents; their contributions accumulate.
    """
    overrides = system.stutter_overrides or {}
    st_lus = float(overrides.get("lus", default_lus))
    st_slus = float(overrides.get("slus", default_slus))
    pf = float(overrides.get("plus_factor", plus_factor))
    motifs = [m for m in system.motif.split(",") if m]
    if not motifs:
        return {}

    expected: dict[str, float] = {}
    for parent in parents:
        cov = parent.n_reads
        if cov <= 0 or not parent.consensus:
            continue
        for variant_seq, step, source in virtual_stutters(parent.consensus, motifs):
            base = st_lus if source == "LUS" else st_slus
            sign_factor = pf if _is_plus_stutter(parent.consensus, variant_seq) else 1.0
            expected[variant_seq] = expected.get(variant_seq, 0.0) + (base**step) * cov * sign_factor
    return expected


def _is_plus_stutter(parent_seq: str, variant_seq: str) -> bool:
    """A ``+`` stutter makes the variant longer than the parent."""
    return len(variant_seq) > len(parent_seq)
