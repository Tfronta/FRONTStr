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
- ``INEXACT_ALLELE`` — an allele whose sequence is a reconstruction rather than
  an observation, so anything sequence-level (ISFG string, iso-allele,
  microvariant) is provisional. Only LongTR ever produced these; with LongTR
  unwired the condition cannot currently arise, and the check remains so the
  flag works if a reconstructing caller is wired in again.
- ``ALLELE_IMBALANCE`` — a heterozygote whose two alleles carry very different
  coverage. The call itself stands; what it warns about is the same locus, or a
  similar one, dropping its minor allele elsewhere in the run.
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
    #: Measured against :attr:`~frontstr.interp.models.MarkerResult.called_reads`
    #: — the reads supporting the genotype — not every read spanning the window.
    #: Reads that clustered into neither allele are not draws from the pair the
    #: derivation models. Using the spanning total instead hid a real dropout:
    #: HG00263 D18S51 is called on 11 reads out of 33 spanning, misses the
    #: second allele Illumina sees, and raised no flag.
    #:
    #: **Derivation, re-run against ONT data (2026-07).** The risk is a true
    #: heterozygote whose minor allele falls under ``calling_thresh`` and is
    #: therefore not called, so a homozygote is reported — a false exclusion.
    #: Take the most unbalanced heterozygote the caller still accepts
    #: (``min_phr_for_het`` = 0.4, so the minor allele is 0.4/1.4 = 28.6% of the
    #: pair) and ask how often it lands below ``calling_thresh`` × the spanning
    #: total. Called coverage is a median 0.795 of spanning across 117 called
    #: loci in the reference slices, which converts one to the other.
    #:
    #: The curve is a sawtooth — the read floor is an integer, so the risk jumps
    #: whenever it crosses one — and reading it at a single N is misleading. Read
    #: as "from this coverage upward the risk never again exceeds":
    #:
    #:     floor 12 → 12.2% risk,  3% of loci flagged
    #:     floor 15 → 12.2% risk,  8%
    #:     floor 17 →  9.7% risk, 12%
    #:     floor 20 →  5.7% risk, 25%   ← the knee
    #:     floor 22 →  5.7% risk, 31%
    #:     floor 25 →  4.6% risk, 43%
    #:
    #: 20 is where the curve turns: 17 → 20 halves the risk for 13 points of
    #: flag rate, while 20 → 25 buys one point of risk for eighteen. So the
    #: value survives the re-derivation — but the **claim attached to it did
    #: not**. This docstring previously said the risk at N=20 was ~1.1%; that
    #: figure came from measuring against the spanning total and is wrong under
    #: the corrected model. It is ~5.7%.
    #:
    #: A quarter of loci flagged on ~30x ONT is therefore the expected rate, not
    #: a symptom. Median called coverage in the reference slices is 26 and the
    #: lower quartile is 20, so the floor sits at Q1 by construction.
    low_coverage_reads: int = Field(default=20, ge=0)

    #: Largest :attr:`~frontstr.interp.models.MarkerResult.allele_balance` a
    #: heterozygote may have and still count as balanced. On that scale 0.5 is
    #: perfect and 0.714 is ``min_phr_for_het`` = 0.4, the point below which no
    #: het is called on read counts at all — so this band sits *inside* the
    #: callable range and warns rather than rejects.
    balanced_ab_max: float = Field(default=0.65, gt=0.5, lt=1.0)

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
        _flag_allele_imbalance(result, thr)
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
    if result.called_reads < thr.low_coverage_reads:
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
            f"{result.marker_name} called on {result.called_reads} read(s) supporting "
            f"the genotype (of {result.total_reads} spanning the window), below the "
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


def _flag_allele_imbalance(result: MarkerResult, thr: QcThresholds) -> None:
    """Flag a heterozygote whose two alleles are unevenly covered.

    Deliberately silent when the call was rescued by phasing: there
    ``HP_RESCUED_HET`` already says the read ratio was the problem and that
    haplotype evidence overrode it. Raising both would put two warnings on one
    phenomenon and train a reviewer to skim them.
    """
    ab = result.allele_balance
    if ab is None or ab <= thr.balanced_ab_max:
        return
    if any(a.hp_rescued for a in result.alleles_called):
        return
    major, minor = sorted(result.alleles_called, key=lambda a: -a.n_reads_total)
    _add(
        result,
        FlagCode.ALLELE_IMBALANCE,
        f"{result.marker_name} heterozygote is unevenly covered: allele balance "
        f"{ab:.2f} (above {thr.balanced_ab_max:.2f}), "
        f"{major.number_label} on {major.n_reads_total} reads vs "
        f"{minor.number_label} on {minor.n_reads_total}. Treat the minor allele as "
        "at risk of dropping out elsewhere in this sample.",
    )


def _flag_inexact(result: MarkerResult) -> None:
    inexact = [a for a in result.alleles_called if a.status == AlleleStatus.INEXACT_ALLELE]
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
    # The note describes how the *bracket count* deviates from the kit. When the
    # number came from STRNaming it is the ISFG designation already, so the
    # warning would be actively misleading — it would tell a reviewer not to
    # compare a number that is precisely the one they should compare.
    if all(a.number_method == "strnaming" for a in result.alleles_called):
        return
    _add(
        result,
        FlagCode.CE_NOMENCLATURE_OFFSET,
        f"{result.marker_name}: the reported allele number is a sequence-derived "
        f"repeat count, not the legacy CE-kit designation. {note} "
        "Do not compare this number directly against a CE profile.",
    )
