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

from dataclasses import dataclass, field
from enum import StrEnum

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


@dataclass(slots=True)
class Allele:
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
    expected_stutter: float = 0.0
    status: AlleleStatus = AlleleStatus.PENDING
    longtr_match: bool = False
    longtr_inexact: bool = False
    longtr_bp_diff: int | None = None

    def fraction(self, total_reads: int) -> float:
        if total_reads <= 0:
            return 0.0
        return self.n_reads_total / total_reads


@dataclass(slots=True)
class MarkerResult:
    """All the forensic decisions for one marker, ready for export/report."""

    marker_name: str
    system: System
    alleles: list[Allele]
    alleles_called: list[Allele]
    call_rule: CallRule
    tri_type: TriType
    total_reads: int
    expected_stutter: dict[str, float] = field(default_factory=dict)
    analytical_thresh: float = 0.02
    calling_thresh: float = 0.10
    longtr_result: LongTRResult | None = None
    discordant: bool = False
    warnings: list[str] = field(default_factory=list)
