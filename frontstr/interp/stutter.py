"""Stutter expectation model (LUS / SLUS).

For each parent allele observed in the reads we build a set of **virtual
stutter sequences** by manipulating the longest and second-longest
uninterrupted stretches of the motif (LUS and SLUS). Each virtual stutter
carries an *expected coverage* = ``rate * parent_coverage``.

We emit -1, -2 and +1 stutters; isometric variants (substitutions inside
the LUS) are skipped because LongTR's INEXACT_ALLELE flag is a stronger
forensic signal.

The function returns a mapping ``virtual_sequence -> expected_coverage``.
When a real cluster's consensus matches one of these keys, its expected
stutter is the value at that key (see :mod:`frontstr.interp.classify`).

Where the rates come from
-------------------------

They are measured, not assumed. :mod:`frontstr.panel.stutter_calib` fits
``rate(-1) = max(0, intercept + slope * LUS)`` to real data, with -2 and +1 as
multipliers of that rate.

This replaces the flat toaSTR-inherited constants (10% LUS, 5% SLUS, forward
stutter at half the reverse rate), which are CE / Illumina figures. Measured on
ONT R10 they are wrong in ways that matter: the flat 10% overestimates the -1
step by 2.3x on average while *missing* the dominant effect, which is that the
rate scales with LUS — from 0.012 at LUS 11 to 0.122 at LUS 14, a 10x range
that no single constant can cover. Forward stutter also runs at 0.73 of the
reverse rate rather than 0.5, because on PCR-free ONT this signal is largely
sequencing error inside the repeat array, which is more symmetric than
polymerase slippage.

See ``docs/stutter_calibration.md`` for the measurement and its caveats — in
particular that the shipped model is fitted on PCR-free WGS and therefore has
no PCR slippage component in it.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from frontstr.evidence.cluster import Cluster
from frontstr.motifs import MotifRun, find_motif_runs, top_two_runs
from frontstr.panel.models import System
from frontstr.panel.stutter_calib import DEFAULT_STUTTER_MODEL, StutterModel

__all__ = [
    "DEFAULT_STUTTER_MODEL",
    "MotifRun",
    "StutterModel",
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
    too; :func:`build_expected_stutter` recovers the sign from the length.
    """
    for variant, _run, step, source in _virtual_stutters_with_run(parent_seq, motifs):
        yield variant, abs(step), source


def _virtual_stutters_with_run(
    parent_seq: str, motifs: Iterable[str]
) -> Iterator[tuple[str, MotifRun, int, str]]:
    """As :func:`virtual_stutters`, but keeping the signed step and source run.

    The rate depends on the length of the run that actually slipped, so the
    expectation model needs the run itself, not just which of the two it was.
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
            yield variant, run, delta, source


def build_expected_stutter(
    parents: list[Cluster],
    system: System,
    *,
    model: StutterModel | None = None,
) -> dict[str, float]:
    """Compute expected stutter coverage per virtual stutter sequence.

    The rate is a function of the length of the run that slipped, so a parent
    whose LUS is 14 units gets a far higher expectation than one at 11 — the
    effect a flat per-marker constant cannot express.

    Args:
        parents: Strong cluster candidates (above ``calling_thresh``).
        system: The marker definition. ``system.stutter_overrides`` may carry
            per-locus tuning, all optional:
            ``log_intercept`` / ``log_slope`` (replace the fitted curve),
            ``slus_factor``, and ``step_-1`` / ``step_-2`` / ``step_1``
            (replace a step multiplier). The legacy ``lus`` key is still
            honoured as a flat, LUS-independent rate for labs that validated
            against it.
        model: Calibrated model to use. Defaults to
            :data:`~frontstr.panel.stutter_calib.DEFAULT_STUTTER_MODEL`.

    Returns:
        Mapping ``virtual_sequence -> expected_coverage``. The same key may
        be hit by multiple parents; their contributions accumulate.
    """
    motifs = [m for m in system.motif.split(",") if m]
    if not motifs:
        return {}
    effective = _apply_overrides(model or DEFAULT_STUTTER_MODEL, system)
    flat_rate = _legacy_flat_rate(system)

    expected: dict[str, float] = {}
    for parent in parents:
        cov = parent.n_reads
        if cov <= 0 or not parent.consensus:
            continue
        for variant_seq, run, step, source in _virtual_stutters_with_run(
            parent.consensus, motifs
        ):
            if flat_rate is not None:
                rate = flat_rate ** abs(step) * effective.step_factors.get(str(step), 1.0)
            else:
                rate = effective.rate(run.n_copies, step)
            if source == "SLUS":
                rate *= effective.slus_factor
            if rate <= 0:
                continue
            expected[variant_seq] = expected.get(variant_seq, 0.0) + rate * cov
    return expected


def _apply_overrides(model: StutterModel, system: System) -> StutterModel:
    """Return ``model`` with any per-marker ``stutter_overrides`` folded in."""
    ov = system.stutter_overrides or {}
    if not ov:
        return model
    step_factors = dict(model.step_factors)
    for step in ("-1", "-2", "1"):
        if f"step_{step}" in ov:
            step_factors[step] = float(ov[f"step_{step}"])
    return model.model_copy(
        update={
            "log_intercept": float(ov.get("log_intercept", model.log_intercept)),
            "log_slope": float(ov.get("log_slope", model.log_slope)),
            "slus_factor": float(ov.get("slus_factor", model.slus_factor)),
            "step_factors": step_factors,
        }
    )


def _legacy_flat_rate(system: System) -> float | None:
    """The pre-calibration flat ``lus`` override, if a panel still sets one."""
    ov = system.stutter_overrides or {}
    return float(ov["lus"]) if "lus" in ov else None
