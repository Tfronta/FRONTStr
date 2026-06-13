"""Interpretation-layer data model.

The Interpretation layer is forensically opinionated. It takes:

- :class:`frontstr.evidence.cluster.Cluster` — what the reads actually say
- Optionally :class:`frontstr.caller.vcf.LongTRResult` — what LongTR called

…and emits :class:`Allele` (per cluster, with classification) + a
:class:`MarkerResult` summarising the called genotype.

These dataclasses are the **forensic source of truth** consumed by exports
and the HTML report. They never reference cyvcf2 or pysam.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field, computed_field

from frontstr.caller.vcf import LongTRResult
from frontstr.panel.models import System


class AlleleStatus(StrEnum):
    """How one cluster was classified by :mod:`frontstr.interp.classify`."""

    PENDING = "pending"
    NO_DATA = "no_data"
    DELETION = "deletion"
    NOISE = "noise"
    STUTTER = "stutter"
    ARTEFACT = "artefact"
    INEXACT_ALLELE = "inexact_allele"
    ALLELE = "allele"


class TriType(StrEnum):
    """Triallelic / mixture pattern for a marker call.

    Empty string is the "regular" homo/het case where triallelism does not
    apply. ``mixture_suspected`` is the default for >2 alleles at a locus
    without ``allow_triallelic=True``.
    """

    NONE = ""
    TYPE_I_UNBALANCED = "tri_I_unbalanced"
    TYPE_II_BALANCED = "tri_II_balanced"
    MIXTURE_SUSPECTED = "mixture_suspected"


class CallRule(StrEnum):
    """Top-level outcome of :func:`frontstr.interp.triallelic.call_profile`."""

    NO_DATA = "no_data"
    HOMOZYGOUS = "homozygous"
    HETEROZYGOUS = "heterozygous"
    TRIALLELIC_TYPE_I = "triallelic_type_I"
    TRIALLELIC_TYPE_II = "triallelic_type_II"
    TRIALLELIC_REVIEW = "triallelic_review"
    TWO_CALLED_THREE_PRESENT = "two_called_three_present_review"


class FlagSeverity(StrEnum):
    """Severity of a structured :class:`Flag`."""

    INFO = "info"
    WARN = "warn"
    ERROR = "error"


class FlagCode(StrEnum):
    """Enumerated, machine-readable review/QC conditions.

    Replaces free-text warning strings so flags can be filtered, aggregated
    (e.g. "all dropouts across a 150-sample batch") and audited. New codes are
    additive; never repurpose an existing string value (it is part of the
    serialized contract).
    """

    LOW_COVERAGE = "low_coverage"
    DROPOUT = "dropout"
    STRAND_BIAS = "strand_bias"
    TRIALLELIC = "triallelic"
    MIXTURE_SUSPECTED = "mixture_suspected"
    LONGTR_DISCORDANT = "longtr_discordant"
    INEXACT_ALLELE = "inexact_allele"
    ISOALLELE = "isoallele"
    CE_NOMENCLATURE_OFFSET = "ce_nomenclature_offset"


_DEFAULT_SEVERITY: dict[FlagCode, FlagSeverity] = {
    FlagCode.LOW_COVERAGE: FlagSeverity.WARN,
    FlagCode.DROPOUT: FlagSeverity.WARN,
    FlagCode.STRAND_BIAS: FlagSeverity.WARN,
    FlagCode.TRIALLELIC: FlagSeverity.INFO,
    FlagCode.MIXTURE_SUSPECTED: FlagSeverity.WARN,
    FlagCode.LONGTR_DISCORDANT: FlagSeverity.WARN,
    FlagCode.INEXACT_ALLELE: FlagSeverity.INFO,
    FlagCode.ISOALLELE: FlagSeverity.INFO,
    FlagCode.CE_NOMENCLATURE_OFFSET: FlagSeverity.INFO,
}


class Flag(BaseModel):
    """One structured, auditable review/QC condition on an allele or marker."""

    code: FlagCode
    severity: FlagSeverity
    message: str

    @classmethod
    def of(
        cls, code: FlagCode, message: str, severity: FlagSeverity | None = None
    ) -> Flag:
        """Build a flag, defaulting severity from the code when not given."""
        return cls(
            code=code,
            severity=severity or _DEFAULT_SEVERITY.get(code, FlagSeverity.INFO),
            message=message,
        )


class IsoAllele(BaseModel):
    """ISFG iso-allele annotation: same allele number, different sequence.

    Folds the former flat ``catalog_*`` fields into one object. ``suffix`` is
    the clean ISFG iso-allele letter ("a"/"b"/…) with no trailing marker;
    ``match_type`` ("exact" | "approx" | "none") replaces the old ``*`` overload.
    ``is_isoallele`` is set at marker level when this allele shares its number
    with a sibling of different sequence (or the catalog resolved a variant).
    """

    suffix: str | None = None
    match_type: str = "none"
    distance: int | None = None
    source: str = ""
    is_isoallele: bool = False


class Allele(BaseModel):
    """One forensic allele candidate. There is exactly one per evidence cluster."""

    cluster_index: int
    consensus: str
    length_bp: int
    n_reads_total: int
    n_reads_hp1: int
    n_reads_hp2: int
    n_reads_hp_none: int
    n_forward: int
    n_reverse: int
    mean_qual: float
    ce: float | None
    isfg: str
    bp_diff: int
    is_deletion: bool
    #: Primary numeric allele for reports: CE when period is defined, else
    #: LongTR-style offset from REF anchor (see :mod:`frontstr.interp.allele_numeric`).
    allele_numeric: float | None = None
    #: ``period_ce`` | ``reference_offset`` | ``delta_only`` | ``deletion`` | ``unavailable``
    allele_numeric_source: str = ""
    expected_stutter: float = 0.0
    status: AlleleStatus = AlleleStatus.PENDING
    longtr_match: bool = False
    longtr_inexact: bool = False
    longtr_bp_diff: int | None = None
    #: ISFG iso-allele annotation (catalog match + same-number sibling). Folds
    #: the former flat ``catalog_suffix``/``catalog_distance``/``catalog_source``.
    iso: IsoAllele = Field(default_factory=IsoAllele)
    #: Structured, auditable per-allele conditions (strand bias, inexact, …).
    flags: list[Flag] = Field(default_factory=list)

    def _number_and_method(self) -> tuple[float | None, str]:
        """Resolve the single canonical allele number + how it was derived.

        Precedence (decided once here, not in the report layer):

        - ``period_ce`` / ``reference_offset`` — calibrated absolute CE number.
        - ``bracket_count`` — sequence-derived repeat count (``ce_from_brackets``)
          for compound markers; the cross-comparable default when there is no
          curated ``reference_ce``.
        - ``delta`` — last-resort relative offset when a compound allele has no
          countable repeat structure.
        - ``bp_sizing`` — raw tandem-repeat span when no allele number exists.
        - ``none`` — deletion / no data.
        """
        src = self.allele_numeric_source
        if self.allele_numeric is not None and src in {"period_ce", "reference_offset"}:
            return float(self.allele_numeric), src
        if src == "delta_only" and self.ce is not None:
            return float(self.ce), "bracket_count"
        if self.allele_numeric is not None and src == "delta_only":
            return float(self.allele_numeric), "delta"
        if self.ce is not None:
            return float(self.ce), "bracket_count"
        if self.length_bp > 0:
            return float(self.length_bp), "bp_sizing"
        return None, "none"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def number(self) -> float | None:
        """Canonical absolute allele number (CE / repeat count); ``None`` if undefined."""
        return self._number_and_method()[0]

    @computed_field  # type: ignore[prop-decorator]
    @property
    def number_method(self) -> str:
        """How :attr:`number` was derived (period_ce | bracket_count | …)."""
        return self._number_and_method()[1]

    def fraction(self, total_reads: int) -> float:
        if total_reads <= 0:
            return 0.0
        return self.n_reads_total / total_reads


class MarkerResult(BaseModel):
    """All the forensic decisions for one marker, ready for export/report."""

    # LongTRResult is a stdlib dataclass; treat it as an opaque nested object
    # rather than letting Pydantic re-validate/copy the caller's internals.
    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: Version of this serialized record's shape. Bump on breaking changes to
    #: the canonical schema so downstream consumers can branch defensively.
    schema_version: str = "1.0"

    marker_name: str
    system: System
    alleles: list[Allele]
    alleles_called: list[Allele]
    call_rule: CallRule
    tri_type: TriType
    total_reads: int
    expected_stutter: dict[str, float] = Field(default_factory=dict)
    analytical_thresh: float = 0.02
    calling_thresh: float = 0.10
    longtr_result: LongTRResult | None = None
    discordant: bool = False
    #: Structured, auditable marker-level conditions. Replaces free-text
    #: warnings so the run is filterable/aggregatable for batch + audit.
    flags: list[Flag] = Field(default_factory=list)
