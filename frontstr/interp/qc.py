"""Run-level QC flags.

:mod:`frontstr.interp.flags` derives the conditions that are intrinsic to a
single finished marker (triallelic, mixture, iso-allele, unpolished consensus).
This module derives the ones that need a *threshold* — a policy decision about
how much coverage is enough, or how skewed a strand ratio has to be before it
stops being chance. Those are laboratory parameters, so they live in
:class:`QcThresholds`, travel into the audit record, and are never hardcoded at
the point of use.

The five codes here were declared in :class:`FlagCode` from the start and never
emitted. A flag that is defined but never raised is worse than no flag at all:
a reviewer scanning a clean report reasonably concludes the condition was
checked and found absent.

Failure modes, in the forensic sense
------------------------------------

- ``DROPOUT`` — the locus produced no call. Reported as a marker-level
  condition, distinct from low coverage: a locus that yields nothing is a
  different problem from one that yields a call on thin evidence.
- ``LOW_COVERAGE`` — a call was made, but on few enough reads that the minor
  allele of a heterozygote could plausibly have been missed. This is the
  allele-dropout risk, and it is why the flag fires on *called* loci.
- ``STRAND_BIAS`` — an allele's reads come overwhelmingly from one strand.
  In a tandem repeat on ONT this can mean a strand-specific basecalling
  artefact rather than a real allele.
- ``INEXACT_ALLELE`` — the caller could not reconstruct the allele exactly, so
  the sequence is a reconstruction. Anything sequence-level (ISFG string,
  iso-allele, microvariant) is provisional.
- ``CE_NOMENCLATURE_OFFSET`` — the reported allele number for this marker is
  known not to equal the legacy CE-kit designation. Curated per marker in the
  panel, because which markers diverge is a property of the kit convention a
  laboratory compares against, not of this code.
"""

from __future__ import annotations

import math

from pydantic import BaseModel, Field

from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    Flag,
    FlagCode,
    MarkerResult,
)


class QcThresholds(BaseModel):
    """Laboratory QC policy. Serialized into the audit record verbatim."""

    #: A called locus below this many reads risks allele dropout.
    #:
    #: The default is derived rather than chosen. Take the most unbalanced
    #: heterozygote the caller still accepts (``min_phr_for_het`` = 0.4, so the
    #: minor allele is 0.4/1.4 = 28.6% of reads) and ask how often it lands
    #: below ``calling_thresh``. Binomial, by coverage:
    #:
    #:     N=10 → 3.5%   N=12 → 10.2%   N=15 → 4.5%
    #:     N=20 → 1.1%   N=25 →  1.3%   N=30 → 0.3%
    #:
    #: 20 is where that risk settles at ~1% and stays there, so it is the point
    #: below which a homozygous call genuinely deserves a second look. (The
    #: non-monotonicity below 20 is the discreteness of the read-count floor,
    #: not noise in the estimate.)
    low_coverage_reads: int = Field(default=20, ge=0)

    #: Two-sided p-value below which a strand ratio is called biased. Kept
    #: strict: at ONT panel coverages a 5% cutoff fires on ordinary sampling
    #: noise often enough to train reviewers to ignore the flag.
    strand_bias_p: float = Field(default=0.01, gt=0.0, le=1.0)

    #: Alleles with fewer reads than this are not strand-tested at all. Below
    #: it the exact binomial cannot reach ``strand_bias_p`` even for a perfect
    #: 0/n split, so testing would only ever produce false reassurance.
    strand_bias_min_reads: int = Field(default=10, ge=2)


def two_sided_binomial_p(successes: int, n: int) -> float:
    """Exact two-sided binomial p-value against a fair coin.

    Written out rather than pulled from SciPy: this is the only statistical
    test in the codebase, read counts are small, and an exact computation is
    both auditable and one fewer compiled dependency.

    Exact because the null is symmetric — ``p = 0.5`` means the two tails have
    equal mass, so doubling the smaller tail is the exact two-sided value
    rather than an approximation.
    """
    if n <= 0:
        return 1.0
    k = min(successes, n - successes)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2**n)
    return float(min(1.0, 2.0 * tail))


def derive_run_qc_flags(
    results: list[MarkerResult], thresholds: QcThresholds | None = None
) -> QcThresholds:
    """Append threshold-dependent QC flags to every marker, in place.

    Idempotent per code, like :func:`frontstr.interp.flags.derive_marker_flags`:
    a code already present on a marker is not added twice, so this may run
    after other flag producers.

    Returns:
        The thresholds actually applied, for the audit record.
    """
    thr = thresholds or QcThresholds()
    for result in results:
        _flag_coverage(result, thr)
        _flag_strand_bias(result, thr)
        _flag_inexact(result)
        _flag_kit_nomenclature(result)
    return thr


def _add(result: MarkerResult, code: FlagCode, message: str) -> None:
    if any(f.code == code for f in result.flags):
        return
    result.flags.append(Flag.of(code, message))


def _flag_coverage(result: MarkerResult, thr: QcThresholds) -> None:
    if result.call_rule == CallRule.NO_DATA:
        _add(
            result,
            FlagCode.DROPOUT,
            f"No allele called at {result.marker_name} "
            f"({result.total_reads} read(s) at the locus).",
        )
        return
    if result.total_reads < thr.low_coverage_reads:
        # A haploid locus has no second allele to drop out, so the usual
        # dropout wording would be misleading; the call is still thin.
        risk = (
            "the call rests on thin evidence"
            if result.system.category == "y_chromosomal"
            else "a minor allele could have been missed, so treat a homozygous "
            "call here as provisional"
        )
        _add(
            result,
            FlagCode.LOW_COVERAGE,
            f"{result.marker_name} called on {result.total_reads} reads, below the "
            f"{thr.low_coverage_reads}-read floor; {risk}.",
        )


def _flag_strand_bias(result: MarkerResult, thr: QcThresholds) -> None:
    biased: list[tuple[Allele, float]] = []
    for a in result.alleles_called:
        n = a.n_forward + a.n_reverse
        if n < thr.strand_bias_min_reads:
            continue
        p = two_sided_binomial_p(a.n_forward, n)
        if p < thr.strand_bias_p:
            biased.append((a, p))
    if not biased:
        return
    detail = ", ".join(
        f"{a.number_label or f'cluster {a.cluster_index}'} "
        f"{a.n_forward}+/{a.n_reverse}- (p={p:.2g})"
        for a, p in biased
    )
    _add(
        result,
        FlagCode.STRAND_BIAS,
        f"Strand-skewed allele(s) at {result.marker_name}: {detail}. "
        "In a tandem repeat this can be a strand-specific basecalling artefact "
        "rather than a real allele.",
    )


def _flag_inexact(result: MarkerResult) -> None:
    inexact = [
        a
        for a in result.alleles_called
        if a.status == AlleleStatus.INEXACT_ALLELE or a.longtr_inexact
    ]
    if not inexact:
        return
    _add(
        result,
        FlagCode.INEXACT_ALLELE,
        f"{len(inexact)} called allele(s) at {result.marker_name} are caller "
        "reconstructions rather than observed sequences; the ISFG string, "
        "iso-allele match and any microvariant are provisional.",
    )


def _flag_kit_nomenclature(result: MarkerResult) -> None:
    """Warn when this marker's number is knowingly not the kit designation.

    Curated in the panel rather than detected here: which markers diverge
    depends on the kit convention a laboratory compares against.
    """
    note = result.system.kit_nomenclature_note
    if not note or not result.alleles_called:
        return
    _add(
        result,
        FlagCode.CE_NOMENCLATURE_OFFSET,
        f"{result.marker_name}: the reported allele number is a sequence-derived "
        f"repeat count, not the legacy CE-kit designation. {note} "
        "Do not compare this number directly against a CE profile.",
    )
