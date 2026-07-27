"""Tests for :mod:`frontstr.panel.stutter_calib`.

The measurement is the deliverable here, so the tests are mostly about the ways
a stutter measurement can be quietly wrong: conditioning on stutter being
present, halving a peak that got split across clusters, measuring at a position
that is really an allele, and extrapolating a fit past its data.
"""

from __future__ import annotations

import math

import pytest

from frontstr.panel.models import System
from frontstr.panel.stutter_calib import (
    DEFAULT_STUTTER_MODEL,
    StutterModel,
    StutterObservation,
    dump_stutter_model,
    fit_stutter_model,
    load_stutter_model,
    lus_units,
    observe_marker,
)

FLANK_L = "GCTTCCGAGTGCAGGTCACAGGGAACACAGACTCCATGGTG"
FLANK_R = "AGGGAAATAAGGGAGGAACAGGCCTTTGGGAATCACCCCAG"


def _system(name: str = "TH01", **kw: object) -> System:
    base: dict[str, object] = {
        "name": name,
        "chromosome": "chr11",
        "ref_start": 1,
        "ref_end": 400,
        "motif": "AATG",
        "period": 4,
    }
    base.update(kw)
    return System(**base)  # type: ignore[arg-type]


def _seq(units: int) -> str:
    return FLANK_L + "AATG" * units + FLANK_R


# ---------------------------------------------------------------------------
# lus_units
# ---------------------------------------------------------------------------


def test_lus_counts_the_longest_run() -> None:
    assert lus_units("AATG" * 9, ["AATG"]) == 9


def test_lus_on_minus_strand_requires_reverse_complement() -> None:
    """Skipping the RC silently returns 0 for every minus-strand marker."""
    from frontstr.motifs import reverse_complement

    rc = reverse_complement("AATG" * 9)
    assert lus_units(rc, ["AATG"], "+") < 9, "sanity: motif is absent as written"
    assert lus_units(rc, ["AATG"], "-") == 9


# ---------------------------------------------------------------------------
# observe_marker
# ---------------------------------------------------------------------------


def test_zero_stutter_positions_are_recorded() -> None:
    """The bias that inflated the first analysis: only counting observed stutter.

    A parent with no stutter at all must still produce observations, otherwise
    the estimate is conditioned on stutter being present.
    """
    out = observe_marker(
        sample="S",
        system=_system(),
        alleles=[(_seq(9), 30)],
        called=[_seq(9)],
    )
    assert {o.step for o in out} == {-1, -2, 1}
    assert all(o.stutter_reads == 0 for o in out)


def test_reads_at_a_position_are_summed_across_clusters() -> None:
    """A stutter peak split into two clusters must not measure as half a peak."""
    out = observe_marker(
        sample="S",
        system=_system(),
        alleles=[(_seq(9), 30), (_seq(8), 2), (_seq(8) + "", 1)],
        called=[_seq(9)],
    )
    minus_one = next(o for o in out if o.step == -1)
    assert minus_one.stutter_reads == 3
    assert minus_one.ratio == pytest.approx(3 / 30)


def test_close_heterozygote_is_excluded() -> None:
    """Alleles 1 unit apart: the -1 position of one IS the other allele."""
    out = observe_marker(
        sample="S",
        system=_system(),
        alleles=[(_seq(9), 30), (_seq(8), 25)],
        called=[_seq(9), _seq(8)],
    )
    assert out == []


def test_well_separated_heterozygote_is_used() -> None:
    out = observe_marker(
        sample="S",
        system=_system(),
        alleles=[(_seq(14), 30), (_seq(8), 25)],
        called=[_seq(14), _seq(8)],
    )
    assert {o.lus for o in out} == {14, 8}


def test_low_coverage_parent_is_skipped() -> None:
    """A ratio from 4 reads carries no information."""
    out = observe_marker(
        sample="S",
        system=_system(),
        alleles=[(_seq(9), 4)],
        called=[_seq(9)],
    )
    assert out == []


def test_parent_reads_come_from_the_position_not_the_cluster() -> None:
    """Coverage at the parent position also sums across split clusters."""
    out = observe_marker(
        sample="S",
        system=_system(),
        alleles=[(_seq(9), 6), (_seq(9), 6), (_seq(8), 3)],
        called=[_seq(9)],
    )
    minus_one = next(o for o in out if o.step == -1)
    assert minus_one.parent_reads == 12
    assert minus_one.stutter_reads == 3


# ---------------------------------------------------------------------------
# fit_stutter_model
# ---------------------------------------------------------------------------


def _obs(lus: int, ratio: float, n: int = 8, parent: int = 10_000) -> list[StutterObservation]:
    """Synthetic observations at a known ratio.

    ``parent`` is large so that rounding reads to integers does not quantise
    the ratio — at 100 reads a true rate of 0.0067 rounds to 0.01, which is
    enough to visibly bend a fit.
    """
    return [
        StutterObservation(
            sample=f"s{i}",
            marker="M",
            step=-1,
            lus=lus,
            parent_reads=parent,
            stutter_reads=round(ratio * parent),
        )
        for i in range(n)
    ]


def test_fit_recovers_a_known_log_linear_trend() -> None:
    obs: list[StutterObservation] = []
    for lus in (10, 11, 12, 13, 14):
        obs += _obs(lus, math.exp(-12.0 + 0.7 * lus))
    model = fit_stutter_model(obs)
    assert model.log_slope == pytest.approx(0.7, abs=0.05)
    assert model.r_squared is not None and model.r_squared > 0.95


def test_fit_range_comes_from_the_bins_that_survived() -> None:
    obs = _obs(11, 0.012) + _obs(13, 0.060) + _obs(4, 0.5, n=1)
    model = fit_stutter_model(obs)
    assert (model.lus_min, model.lus_max) == (11, 13)


def test_thin_bins_are_excluded_from_the_fit() -> None:
    """A single noisy parent at an extreme LUS has huge regression leverage."""
    clean = _obs(11, 0.012) + _obs(13, 0.060)
    with_outlier = clean + _obs(4, 0.9, n=1)
    assert fit_stutter_model(with_outlier).log_slope == pytest.approx(
        fit_stutter_model(clean).log_slope
    )


def test_fit_refuses_rather_than_guessing_when_data_is_too_thin() -> None:
    with pytest.raises(ValueError, match="at least 2 usable LUS bins"):
        fit_stutter_model(_obs(11, 0.012))


def test_fit_needs_minus_one_observations() -> None:
    with pytest.raises(ValueError, match="no -1 step"):
        fit_stutter_model([])


def test_step_factors_are_ratios_of_the_minus_one_rate() -> None:
    obs = _obs(11, 0.02) + _obs(13, 0.06)
    obs += [
        StutterObservation(
            sample="s", marker="M", step=-2, lus=11, parent_reads=100, stutter_reads=1
        ),
        StutterObservation(
            sample="s", marker="M", step=1, lus=11, parent_reads=100, stutter_reads=2
        ),
    ]
    model = fit_stutter_model(obs)
    pooled_minus1 = sum(o.stutter_reads for o in obs if o.step == -1) / sum(
        o.parent_reads for o in obs if o.step == -1
    )
    assert model.step_factors["-2"] == pytest.approx(0.01 / pooled_minus1, abs=1e-3)


# ---------------------------------------------------------------------------
# StutterModel.rate
# ---------------------------------------------------------------------------


def test_rate_is_clamped_outside_the_calibrated_range() -> None:
    m = DEFAULT_STUTTER_MODEL
    assert m.rate(2, -1) == pytest.approx(m.rate(m.lus_min, -1))
    assert m.rate(99, -1) == pytest.approx(m.rate(m.lus_max, -1))


def test_rate_is_never_negative_and_grows_with_lus() -> None:
    m = DEFAULT_STUTTER_MODEL
    rates = [m.rate(lus, -1) for lus in range(m.lus_min, m.lus_max + 1)]
    assert all(r > 0 for r in rates)
    assert rates == sorted(rates)


def test_shipped_model_matches_the_documented_measurements() -> None:
    """Guards docs/stutter_calibration.md against silent drift."""
    m = DEFAULT_STUTTER_MODEL
    assert m.protocol == "wgs_pcr_free", "the PCR-free caveat must survive"
    for lus, observed in ((10, 0.0100), (12, 0.0346), (14, 0.1222)):
        assert m.rate(lus, -1) == pytest.approx(observed, abs=0.01)


def test_unknown_step_has_no_expectation() -> None:
    assert DEFAULT_STUTTER_MODEL.rate(12, -3) == 0.0


# ---------------------------------------------------------------------------
# Round-trip
# ---------------------------------------------------------------------------


def test_model_round_trips_through_json(tmp_path) -> None:
    path = dump_stutter_model(DEFAULT_STUTTER_MODEL, tmp_path / "m.json")
    loaded = load_stutter_model(path)
    assert loaded == DEFAULT_STUTTER_MODEL
    assert isinstance(loaded, StutterModel)
