"""Interpretation-layer data model.

The Interpretation layer is forensically opinionated. It takes:

- :class:`frontstr.evidence.cluster.Cluster` — what the reads actually say

…and emits :class:`Allele` (per cluster, with classification) + a
:class:`MarkerResult` summarising the called genotype.

These dataclasses are the **forensic source of truth** consumed by exports
and the HTML report. They never reference cyvcf2 or pysam.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, Field, computed_field

from frontstr.evidence.consensus import ConsensusMethod
from frontstr.panel.models import System


class AlleleStatus(StrEnum):
    """How one cluster was classified by :mod:`frontstr.interp.classify`."""

    PENDING = "pending"
    NO_DATA = "no_data"
    DELETION = "deletion"
    NOISE = "noise"
    STUTTER = "stutter"
    ARTEFACT = "artefact"
    #: Same allele as a stronger cluster on the same haplotype, split apart by
    #: sequencing error. See :mod:`frontstr.interp.haplotype`.
    HP_PHANTOM = "hp_phantom"
    #: Retired with LongTR: it marked clusters the caller could only
    #: reconstruct. Nothing sets it now; kept so historical records parse.
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
    #: Retired: LongTR is no longer wired into the pipeline (see
    #: :mod:`frontstr.caller`). Kept so historical records still parse; nothing
    #: produces it.
    LONGTR_DISCORDANT = "longtr_discordant"
    #: Retired with LongTR — see :class:`AlleleStatus.INEXACT_ALLELE`.
    INEXACT_ALLELE = "inexact_allele"
    ISOALLELE = "isoallele"
    CE_NOMENCLATURE_OFFSET = "ce_nomenclature_offset"
    CONSENSUS_FALLBACK = "consensus_fallback"
    HP_PHANTOM_COLLAPSED = "hp_phantom_collapsed"
    HP_RESCUED_HET = "hp_rescued_het"
    PHASE_BLOCK_SPLIT = "phase_block_split"
    ALLELE_IMBALANCE = "allele_imbalance"
    NON_DEFAULT_THRESHOLD = "non_default_threshold"


_DEFAULT_SEVERITY: dict[FlagCode, FlagSeverity] = {
    FlagCode.LOW_COVERAGE: FlagSeverity.WARN,
    FlagCode.DROPOUT: FlagSeverity.WARN,
    FlagCode.STRAND_BIAS: FlagSeverity.WARN,
    FlagCode.TRIALLELIC: FlagSeverity.INFO,
    FlagCode.MIXTURE_SUSPECTED: FlagSeverity.WARN,
    FlagCode.LONGTR_DISCORDANT: FlagSeverity.WARN,
    FlagCode.INEXACT_ALLELE: FlagSeverity.INFO,
    FlagCode.ISOALLELE: FlagSeverity.INFO,
    # WARN, not INFO: comparing this number against a CE profile without
    # knowing it is a bracket count can produce a false exclusion.
    FlagCode.CE_NOMENCLATURE_OFFSET: FlagSeverity.WARN,
    # WARN, not INFO: an unpolished consensus carries single-read errors into
    # the ISFG string and the iso-allele catalog match.
    FlagCode.CONSENSUS_FALLBACK: FlagSeverity.WARN,
    # INFO: the collapse is the *correct* call, but it must stay visible so a
    # reviewer can see that a candidate was suppressed and why.
    FlagCode.HP_PHANTOM_COLLAPSED: FlagSeverity.INFO,
    # INFO, same reasoning as above: the heterozygous call is the correct one,
    # but a reviewer must be able to see that it rests on phasing rather than on
    # peak balance, because the read ratio alone would read as homozygous.
    FlagCode.HP_RESCUED_HET: FlagSeverity.INFO,
    # WARN: haplotype evidence is weaker than it looks at this locus. Every
    # haplotype rule declines rather than guessing, so the call is safe — but a
    # reviewer comparing HP counts by eye would draw a conclusion the caller
    # deliberately refused to draw.
    FlagCode.PHASE_BLOCK_SPLIT: FlagSeverity.WARN,
    # WARN: an imbalanced heterozygote is the visible half of a locus that may
    # be dropping the other allele elsewhere in the run — degradation, a primer
    # variant under a flank, or simply thin coverage.
    FlagCode.ALLELE_IMBALANCE: FlagSeverity.WARN,
    # WARN: the run overrode a threshold whose default was derived from measured
    # data. The call may be perfectly good — the point is that it is not
    # comparable with a default run, and six months later nobody will remember.
    FlagCode.NON_DEFAULT_THRESHOLD: FlagSeverity.WARN,
}


class Flag(BaseModel):
    """One structured, auditable review/QC condition on an allele or marker."""

    code: FlagCode
    severity: FlagSeverity
    message: str

    @classmethod
    def of(cls, code: FlagCode, message: str, severity: FlagSeverity | None = None) -> Flag:
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


def _trim_number(value: float) -> str:
    """Render an allele number without trailing zeros (``9.3``, ``14``)."""
    x = round(float(value), 4)
    if abs(x - int(x)) < 1e-9:
        return str(round(x))
    return f"{x:.10f}".rstrip("0").rstrip(".")


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
    #: STRNaming's full allele name over the marker's standard reporting range,
    #: e.g. ``CE29_TCTA[4]TCTG[6]...``. Empty when STRNaming defines no range for
    #: the marker or the range could not be located in the consensus.
    strnaming_name: str = ""
    #: CE parsed from :attr:`strnaming_name`. This is the canonical allele number
    #: whenever it is present — see :meth:`_number_and_method`.
    strnaming_ce: float | None = None
    #: Why STRNaming did or did not name this allele
    #: (:class:`frontstr.interp.naming.NameStatus`). Kept for audit: a reviewer
    #: must be able to see that a number came from the legacy path and why.
    strnaming_status: str = ""
    bp_diff: int
    is_deletion: bool
    #: How ``consensus`` was derived (``poa_spoa`` | ``poa_abpoa`` | ``single``
    #: | ``mode`` | ``empty``). ``mode`` means no POA backend was available and
    #: the sequence is a single unpolished read — see ``CONSENSUS_FALLBACK``.
    consensus_method: str = ConsensusMethod.SINGLE.value
    #: Primary numeric allele for reports: CE when period is defined, else an
    #: offset from the REF anchor (see :mod:`frontstr.interp.allele_numeric`).
    allele_numeric: float | None = None
    #: ``period_ce`` | ``reference_offset`` | ``delta_only`` | ``deletion`` | ``unavailable``
    allele_numeric_source: str = ""
    expected_stutter: float = 0.0
    #: Reads recovered from same-haplotype phantom clusters of this allele
    #: (see :mod:`frontstr.interp.haplotype`). Reported so coverage is not
    #: understated; ``n_reads_total`` deliberately stays as observed.
    n_reads_absorbed: int = 0
    #: Phase block (BAM ``PS`` tag) this allele's haplotype-tagged reads came
    #: from, or ``None`` when they carry no ``PS`` or span more than one.
    phase_set: int | None = None
    #: How many distinct phase blocks those reads came from. Anything above 1
    #: means the ``HP`` labels within this cluster are not comparable, so
    #: :func:`frontstr.interp.haplotype.dominant_hp` refuses to assign it a
    #: haplotype at all.
    n_phase_sets: int = 0
    #: Set when this allele was called only because phasing put it on the
    #: opposite haplotype from the major allele — the peak-height ratio alone
    #: would have collapsed the locus to homozygous. Raises ``HP_RESCUED_HET``.
    hp_rescued: bool = False
    status: AlleleStatus = AlleleStatus.PENDING
    #: ISFG iso-allele annotation (catalog match + same-number sibling). Folds
    #: the former flat ``catalog_suffix``/``catalog_distance``/``catalog_source``.
    iso: IsoAllele = Field(default_factory=IsoAllele)
    #: Structured, auditable per-allele conditions (strand bias, inexact, …).
    flags: list[Flag] = Field(default_factory=list)

    def _number_and_method(self) -> tuple[float | None, str]:
        """Resolve the single canonical allele number + how it was derived.

        Precedence (decided once here, not in the report layer):

        - ``strnaming`` — the ISFG-recommended designation, computed by
          STRNaming over the marker's standard reporting range. Outranks
          everything else because it is the only method that is right on the
          compound markers: the bracket count is off by an allele-structure
          dependent amount at vWA, D21S11, D2S1338, D1S1656 and D2S441, and the
          panel ``corr_value`` is off at DXS7132. See
          :mod:`frontstr.interp.naming`.
        - ``period_ce`` / ``reference_offset`` — calibrated absolute CE number.
        - ``bracket_count`` — sequence-derived repeat count (``ce_from_brackets``)
          for compound markers; the cross-comparable default when there is no
          curated ``reference_ce``.
        - ``delta`` — last-resort relative offset when a compound allele has no
          countable repeat structure.
        - ``bp_sizing`` — raw tandem-repeat span when no allele number exists.
        - ``none`` — deletion / no data.
        """
        if self.strnaming_ce is not None:
            return float(self.strnaming_ce), "strnaming"
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

    @computed_field  # type: ignore[prop-decorator]
    @property
    def number_is_absolute(self) -> bool:
        """True when :attr:`number` is a real, cross-comparable allele number.

        False for the relative and fallback methods (``delta``, ``bp_sizing``,
        ``none``), which must not be read as a kit CE designation.
        """
        return self.number_method in (
            "strnaming",
            "period_ce",
            "reference_offset",
            "bracket_count",
        )

    @computed_field  # type: ignore[prop-decorator]
    @property
    def repeat_label(self) -> str:
        """The single bracketed-repeat string every view renders.

        Same reasoning as :attr:`number_label`, which exists because the CLI and
        the report used to format allele numbers independently and could show
        ``Δ-2`` and ``14`` for one allele. The strings had drifted the same way:
        ``--trace`` printed STRNaming's ``CE9.3_TGAA[6]TGA[1]TGAA[3]`` while the
        HTML, CSV and XLSX printed :attr:`isfg`, which is
        :func:`~frontstr.interp.isfg.compress_isfg` over the **whole panel
        window** — a hundred lowercase flank bases before the brackets even
        start. Two strings for one allele, in a forensic report.

        STRNaming's name wins when present: it is the ISFG DNA Commission's
        prescribed format (Gettings et al. 2024, Recommendation 2), it is scoped
        to the marker's standard reporting range instead of our extraction
        window, and it carries flanking variants explicitly. The legacy string
        remains the fallback for markers STRNaming has no range for (DYS393,
        AMEL) and is still available as :attr:`isfg` for anything that needs the
        raw window view.
        """
        return self.strnaming_name or self.isfg

    @computed_field  # type: ignore[prop-decorator]
    @property
    def repeat_label_source(self) -> str:
        """``strnaming`` | ``bracket_scan`` | ``none`` — how :attr:`repeat_label` was made."""
        if self.strnaming_name:
            return "strnaming"
        return "bracket_scan" if self.isfg else "none"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def number_label(self) -> str:
        """The single human-facing label for this allele.

        Every view — the CLI table, the HTML report, the CSV exports — renders
        this string. Formatting used to be duplicated per view, which let the
        same allele read as ``Δ-2`` in ``frontstr interpret`` and ``14`` in the
        report. In a forensic context two numbers for one allele is not a
        cosmetic problem, so the label is produced once, here, beside the
        number it describes.
        """
        number = self.number
        if number is None:
            # Non-numeric markers (AMEL X / Y) carry their designation in ISFG.
            return self.isfg or ""
        trimmed = _trim_number(number)
        if self.number_method == "delta":
            return f"Δ{trimmed}" if abs(number) > 1e-9 else trimmed
        if self.number_method == "bp_sizing":
            return f"{int(number)} bp TR"
        return trimmed

    def fraction(self, total_reads: int) -> float:
        if total_reads <= 0:
            return 0.0
        return self.n_reads_total / total_reads


class MarkerResult(BaseModel):
    """All the forensic decisions for one marker, ready for export/report."""

    #: Version of this serialized record's shape. Bump on breaking changes to
    #: the canonical schema so downstream consumers can branch defensively.
    #: 2.0 removed the ``longtr_*`` allele fields, ``longtr_result`` and
    #: ``discordant`` when LongTR was unwired, and added the ``strnaming_*``
    #: fields.
    schema_version: str = "2.0"

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
    #: Structured, auditable marker-level conditions. Replaces free-text
    #: warnings so the run is filterable/aggregatable for batch + audit.
    flags: list[Flag] = Field(default_factory=list)

    @computed_field  # type: ignore[prop-decorator]
    @property
    def allele_balance(self) -> float | None:
        """Coverage share of the strongest called allele. ``None`` unless het.

        **Convention: strongest allele over the sum of called alleles**, so the
        value runs from 0.50 (perfectly balanced) to 1.0 (all reads on one
        allele). Stating that matters — with the *strongest* allele on top the
        scale is one-sided, and a band written as if it were symmetric around
        0.5 would have half of itself unreachable.

        This replaces the peak-height ratio in reported output rather than
        joining it. AB and PHR are the same measurement (``AB = 1/(1+PHR)``),
        and shipping both would put two numbers for one quantity in front of a
        reviewer — the mistake that once had the same allele reading ``Δ-2`` in
        the CLI and ``14`` in the report.

        Where the landmarks fall:

            AB 0.500   perfectly balanced
            AB 0.650   edge of the balanced band (``QcThresholds``)
            AB 0.714   ``min_phr_for_het`` = 0.4; below this no het is called
                       on read counts alone
            AB 0.773   HG00113 D2S1338, called het on phasing alone

        Only defined for a two-allele call: a homozygote has nothing to
        balance, and for a triallelic locus a single ratio would be a fiction.
        """
        if len(self.alleles_called) != 2:
            return None
        counts = sorted((a.n_reads_total for a in self.alleles_called), reverse=True)
        total = counts[0] + counts[1]
        if total <= 0:
            return None
        return round(counts[0] / total, 3)
