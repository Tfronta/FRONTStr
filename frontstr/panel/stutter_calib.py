"""Empirical stutter calibration from real data.

Why this exists
---------------

The stutter rates FRONTStr shipped with (10% per LUS step, 5% per SLUS step,
forward stutter at half the reverse rate) are inherited from CE / Illumina
practice via toaSTR. On ONT they are wrong in both directions, and the errors
are not small — see ``docs/stutter_calibration.md``.

This module measures the rates from data instead. It is deliberately a separate
calibration pass rather than something that runs during a case: the model is a
property of the chemistry and the protocol, not of the sample.

Measurement method
------------------

Standard forensic practice, adapted to sequence data:

- Only loci where a stutter position cannot be confused with a real allele are
  used: homozygotes, or heterozygotes whose alleles are at least
  ``MIN_ALLELE_SEPARATION_UNITS`` repeat units apart.
- Reads are grouped by **repeat-core length**, and every read at a given offset
  is summed. A stutter peak is a *position*, not a cluster, so a peak split
  across two clusters must not halve the measured ratio.
- Parents with fewer than ``MIN_PARENT_READS`` reads are skipped — a ratio
  estimated from 4 reads carries no information.
- Positions with **zero** stutter reads are recorded as zeros. Omitting them
  conditions the estimate on stutter being present and inflates every rate
  (measured on the slice set: 0.098 vs the correct 0.044 for the -1 step).

Model form
----------

Rate is **log-linear** in the parent's LUS (longest uninterrupted motif run,
the standard forensic covariate)::

    rate(-1)   = exp(log_intercept + log_slope * clamp(LUS, lus_min, lus_max))
    rate(step) = rate(-1) * step_factor[step]

A plain linear fit was tried first and rejected. The measured rates accelerate
with LUS (0.010, 0.012, 0.035, 0.060, 0.122 at LUS 10-14), so a straight line
undershoots at both ends and — worse — crosses zero inside the range where
stutter is actually observed, which would mean "no stutter model at all" for
short-LUS loci. The log-linear form is convex by construction, cannot go
negative, and treats the rate as a multiplicative hazard, which is what a
per-unit slippage process is.

Outside the fitted LUS range the LUS is **clamped**, not extrapolated: nothing
measured supports the behaviour of the curve at LUS 4 or LUS 30.

The -2 and +1 steps are expressed as multipliers of the -1 rate rather than as
a geometric decay: on ONT the geometric form underestimates -2 by ~2.5x, because
what is being measured is largely sequencing error inside the repeat array
rather than a multi-cycle polymerase slippage process.

Which way to be wrong
---------------------

Over-predicting stutter classifies a real minor allele as stutter, producing a
silent false homozygote — the analyst never sees the allele that was removed.
Under-predicting lets a stutter peak survive as an allele candidate, which
surfaces as a review flag and is additionally caught by ``min_reads_third`` and
haplotype-aware suppression. The second failure mode is far safer, so where the
data is thin the model is left deliberately conservative.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from pydantic import BaseModel, Field

from frontstr.motifs import find_motif_runs, repeat_core_length, reverse_complement
from frontstr.panel.models import Panel, System

#: Alleles closer than this (in repeat units) make the stutter position
#: ambiguous with a real allele, so the locus is excluded from calibration.
MIN_ALLELE_SEPARATION_UNITS = 3

#: A parent with fewer reads than this cannot support a rate estimate.
MIN_PARENT_READS = 8

#: LUS bins with less support than this are excluded from the fit. A bin built
#: from one parent is noise, and a noisy bin at an extreme LUS drags the slope
#: hard because of its regression leverage.
MIN_BIN_OBSERVATIONS = 7

#: Stutter steps modelled, in repeat units.
MODELLED_STEPS = (-1, -2, 1)


@dataclass(frozen=True, slots=True)
class StutterObservation:
    """One (parent, step) measurement at one locus in one sample."""

    sample: str
    marker: str
    step: int
    lus: int
    parent_reads: int
    stutter_reads: int

    @property
    def ratio(self) -> float:
        return self.stutter_reads / self.parent_reads if self.parent_reads else 0.0


class StutterModel(BaseModel):
    """A calibrated stutter expectation model.

    ``protocol`` is load-bearing, not decoration: a model fitted on PCR-free
    WGS has no PCR slippage component in it, so applying it to amplicon
    casework will under-predict stutter. Consumers should compare it against
    the protocol they are running.
    """

    version: str = "1.0"
    #: Free-text provenance: samples, chemistry, basecaller, date.
    source: str = ""
    #: ``wgs_pcr_free`` | ``amplicon`` | ``unknown``.
    protocol: str = "unknown"

    #: rate(-1) = exp(log_intercept + log_slope * clamp(LUS, lus_min, lus_max))
    log_intercept: float
    log_slope: float
    #: LUS range the curve was actually fitted over. Outside it the LUS is
    #: clamped to the nearest end rather than extrapolated — nothing measured
    #: supports the curve's behaviour at LUS 4 or LUS 30.
    lus_min: int = 1
    lus_max: int = 40
    #: Multipliers on rate(-1), keyed by step as a string ("-1", "-2", "1").
    step_factors: dict[str, float] = Field(default_factory=dict)
    #: Extra multiplier applied to stutter originating in the second-longest
    #: run. Defaults to 1.0 because the rate is already a function of the
    #: *slipping run's* length, so a shorter secondary run gets a lower rate
    #: automatically — the old flat ``slus = lus / 2`` was a crude stand-in for
    #: exactly that. Left as a knob for labs that need to tune it; not derived
    #: from data, since the -1 position is the same whichever run produced it
    #: and the two contributions are not separable by this measurement.
    slus_factor: float = 1.0

    n_observations: int = 0
    n_loci: int = 0
    r_squared: float | None = None

    def rate(self, lus: int, step: int) -> float:
        """Expected stutter rate for a parent with ``lus`` units at ``step``.

        ``lus`` is clamped to the calibrated range before the curve is
        evaluated, so a locus outside that range gets the nearest supported
        rate instead of an extrapolation.
        """
        lus_eff = min(max(lus, self.lus_min), self.lus_max)
        base = math.exp(self.log_intercept + self.log_slope * lus_eff)
        return base * self.step_factors.get(str(step), 0.0)

    def describe(self) -> str:
        """One line naming this model, for a run header or a log.

        Carries ``protocol`` because that is the field that invalidates a run:
        this model has no PCR slippage component, so applying it to amplicon
        casework under-predicts stutter. A trace that does not say which model
        decided which candidates were artefacts cannot be reproduced.
        """
        parts = [self.version, self.protocol]
        if self.r_squared is not None:
            parts.append(f"R² {self.r_squared:.3f}")
        if self.n_observations:
            parts.append(f"n={self.n_observations} over {self.n_loci} loci")
        parts.append(f"LUS {self.lus_min}–{self.lus_max}, clamped outside")
        return ", ".join(parts)


#: First-pass model, fitted on the 5 ONT R10/Dorado 1000G slices.
#: See docs/stutter_calibration.md. PCR-free WGS — NOT valid for amplicon.
DEFAULT_STUTTER_MODEL = StutterModel(
    version="2026.07-ont-r10-wgs",
    source=(
        "5 ONT R10.4.1/LSK114 Dorado 1000G WGS slices "
        "(HG00097, HG00113, HG00154, HG00263, GM19038), CODIS 20 + sex panel"
    ),
    protocol="wgs_pcr_free",
    # Reproduce with:
    #   frontstr calibrate-stutter -p examples/panels/codis_20_grch38.yaml \
    #       --protocol wgs_pcr_free --bam tests/data/ont_slices/*.codis.bam
    # Fitted rates vs observed, LUS 10-14:
    #   0.0070/0.0100  0.0145/0.0121  0.0296/0.0346  0.0605/0.0604  0.1237/0.1222
    log_intercept=-12.1125,
    log_slope=0.7159,
    # Only LUS bins with >= MIN_BIN_OBSERVATIONS support survive the fit; on
    # this slice set that is 10-14. Loci outside the range are clamped to the
    # nearest end. Widening this range is the main reason to calibrate on more
    # samples.
    lus_min=10,
    lus_max=14,
    step_factors={"-1": 1.0, "-2": 0.242, "1": 0.726},
    slus_factor=1.0,
    n_observations=76,
    n_loci=52,
    r_squared=0.965,
)


def _canonical(seq: str, strand: str) -> str:
    return reverse_complement(seq) if strand == "-" else seq


def lus_units(sequence: str, motifs: list[str], strand: str = "+") -> int:
    """Longest uninterrupted motif run in ``sequence``, in copies.

    Computed on the canonical strand: a minus-strand marker's motif does not
    appear in the reference-oriented consensus at all, so skipping the
    reverse-complement silently returns 0 for half the panel.
    """
    runs = find_motif_runs(_canonical(sequence, strand), motifs)
    return max((r.n_copies for r in runs), default=0)


def observe_marker(
    *,
    sample: str,
    system: System,
    alleles: list[tuple[str, int]],
    called: list[str],
    min_parent_reads: int = MIN_PARENT_READS,
) -> list[StutterObservation]:
    """Collect stutter observations for one marker in one sample.

    Args:
        sample: Sample identifier, carried into each observation.
        system: Marker definition.
        alleles: ``(consensus, n_reads)`` for **every** cluster at the locus,
            called or not — stutter peaks are among the uncalled ones.
        called: Consensus sequences of the alleles that were actually called.
        min_parent_reads: Skip parents with less support than this.

    Returns:
        One observation per (parent, step), including zero-stutter positions.
        Empty when the locus is unusable for calibration.
    """
    motifs = [m for m in system.motif.split(",") if m]
    if not motifs:
        return []
    period = system.period if system.period > 0 else 4

    reads_at: dict[int, int] = defaultdict(int)
    for consensus, n_reads in alleles:
        core = repeat_core_length(consensus, motifs, strand=system.strand)
        if core is not None:
            reads_at[core] += n_reads

    parent_cores: list[tuple[str, int]] = []
    for consensus in called:
        core = repeat_core_length(consensus, motifs, strand=system.strand)
        if core is not None:
            parent_cores.append((consensus, core))
    if not parent_cores:
        return []

    cores = sorted(core for _, core in parent_cores)
    if len(cores) > 1:
        closest = min(b - a for a, b in pairwise(cores))
        if closest < MIN_ALLELE_SEPARATION_UNITS * period:
            return []  # a stutter position would be indistinguishable from an allele

    out: list[StutterObservation] = []
    for consensus, core in parent_cores:
        parent_reads = reads_at[core]
        if parent_reads < min_parent_reads:
            continue
        lus = lus_units(consensus, motifs, system.strand)
        for step in MODELLED_STEPS:
            out.append(
                StutterObservation(
                    sample=sample,
                    marker=system.name,
                    step=step,
                    lus=lus,
                    parent_reads=parent_reads,
                    stutter_reads=reads_at.get(core + step * period, 0),
                )
            )
    return out


def collect_observations(
    bams: list[Path],
    panel: Panel,
    *,
    min_mapq: int = 20,
    reference_fasta: Path | None = None,
) -> list[StutterObservation]:
    """Run the evidence + interpretation path over ``bams`` and measure stutter."""
    from frontstr.interp.profile import interpret_run

    out: list[StutterObservation] = []
    for bam in bams:
        sample = bam.name.split(".")[0]
        results = interpret_run(
            bam=bam, panel=panel, min_mapq=min_mapq, reference_fasta=reference_fasta
        )
        for r in results:
            if r.system.marker_type != "str" or not r.alleles_called:
                continue
            out.extend(
                observe_marker(
                    sample=sample,
                    system=r.system,
                    alleles=[(a.consensus, a.n_reads_total) for a in r.alleles],
                    called=[a.consensus for a in r.alleles_called],
                )
            )
    return out


def _pooled_ratio(rows: list[StutterObservation]) -> float:
    """Read-weighted ratio: total stutter reads over total parent reads."""
    parents = sum(o.parent_reads for o in rows)
    return sum(o.stutter_reads for o in rows) / parents if parents else 0.0


def fit_stutter_model(
    observations: list[StutterObservation],
    *,
    source: str = "",
    protocol: str = "unknown",
    version: str = "custom",
    min_bin_observations: int = MIN_BIN_OBSERVATIONS,
) -> StutterModel:
    """Fit :class:`StutterModel` to measured observations.

    The LUS trend is fitted on the ``-1`` step by weighted least squares of
    ``log(rate)`` on LUS, each bin weighted by its total parent reads. ``-2``
    and ``+1`` become multipliers on the fitted ``-1`` rate.

    LUS bins with fewer than ``min_bin_observations`` measurements, or with a
    pooled rate of zero, are excluded from the fit: a ratio estimated from one
    parent is noise, and a single noisy bin at an extreme LUS has enormous
    leverage on the slope. The surviving bins define ``lus_min``/``lus_max``,
    outside which the model clamps rather than extrapolates.

    Raises:
        ValueError: If too few LUS bins survive to fit a slope.
    """
    minus_one = [o for o in observations if o.step == -1]
    if not minus_one:
        raise ValueError("no -1 step observations to fit")

    by_lus: dict[int, list[StutterObservation]] = defaultdict(list)
    for o in minus_one:
        by_lus[o.lus].append(o)

    points = [
        (lus, ratio, float(sum(o.parent_reads for o in rows)))
        for lus, rows in sorted(by_lus.items())
        if len(rows) >= min_bin_observations and (ratio := _pooled_ratio(rows)) > 0
    ]
    if len(points) < 2:
        raise ValueError(
            f"need at least 2 usable LUS bins (>= {min_bin_observations} observations "
            f"and a non-zero rate) to fit a slope, got {len(points)}. "
            "Calibrate on more samples."
        )

    log_points = [(x, math.log(y), w) for x, y, w in points]
    total_w = sum(w for *_, w in log_points)
    mean_x = sum(w * x for x, _, w in log_points) / total_w
    mean_y = sum(w * y for _, y, w in log_points) / total_w
    denom = sum(w * (x - mean_x) ** 2 for x, _, w in log_points)
    if denom == 0:
        raise ValueError("all usable bins share one LUS value; cannot fit a slope")
    slope = sum(w * (x - mean_x) * (y - mean_y) for x, y, w in log_points) / denom
    intercept = mean_y - slope * mean_x

    ss_tot = sum(w * (y - mean_y) ** 2 for _, y, w in log_points)
    ss_res = sum(w * (y - (intercept + slope * x)) ** 2 for x, y, w in log_points)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None

    base = _pooled_ratio(minus_one)
    factors: dict[str, float] = {}
    for step in MODELLED_STEPS:
        rows = [o for o in observations if o.step == step]
        factors[str(step)] = (_pooled_ratio(rows) / base) if rows and base else 0.0

    return StutterModel(
        version=version,
        source=source,
        protocol=protocol,
        log_intercept=round(intercept, 4),
        log_slope=round(slope, 4),
        lus_min=min(x for x, *_ in points),
        lus_max=max(x for x, *_ in points),
        step_factors={k: round(v, 3) for k, v in factors.items()},
        n_observations=len(minus_one),
        n_loci=len({(o.sample, o.marker) for o in minus_one}),
        r_squared=round(r2, 3) if r2 is not None else None,
    )


def summarise(observations: list[StutterObservation]) -> dict[str, object]:
    """Human-readable breakdown used by the CLI and by the calibration doc."""
    by_step = {
        str(step): {
            "n": sum(1 for o in observations if o.step == step),
            "pooled_ratio": round(_pooled_ratio([o for o in observations if o.step == step]), 4),
            "n_zero": sum(1 for o in observations if o.step == step and o.stutter_reads == 0),
        }
        for step in MODELLED_STEPS
    }
    by_lus: dict[int, list[StutterObservation]] = defaultdict(list)
    by_marker: dict[str, list[StutterObservation]] = defaultdict(list)
    for o in observations:
        if o.step == -1:
            by_lus[o.lus].append(o)
            by_marker[o.marker].append(o)
    return {
        "by_step": by_step,
        "minus1_by_lus": {
            lus: {"n": len(rows), "pooled_ratio": round(_pooled_ratio(rows), 4)}
            for lus, rows in sorted(by_lus.items())
        },
        "minus1_by_marker": {
            m: {"n": len(rows), "pooled_ratio": round(_pooled_ratio(rows), 4)}
            for m, rows in sorted(by_marker.items())
        },
    }


def load_stutter_model(path: Path) -> StutterModel:
    """Load a calibrated model from JSON."""
    return StutterModel.model_validate_json(path.read_text(encoding="utf-8"))


def dump_stutter_model(model: StutterModel, path: Path) -> Path:
    """Write a calibrated model to JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(model.model_dump(mode="json"), indent=2) + "\n", encoding="utf-8")
    return path
