"""Pydantic models for panels and markers (Systems).

A ``System`` is a single STR locus; a ``Panel`` groups multiple systems.

Phase 1.2 of ROADMAP.md. See plan-longtr-improved.md §3.4 and §6.
"""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator


class System(BaseModel):
    """A single forensic STR marker (a.k.a. locus, system)."""

    name: str
    codis_name: str | None = None
    chromosome: str
    ref_start: int = Field(ge=1)
    ref_end: int = Field(ge=1)
    motif: str
    period: int = Field(ge=-1)
    corr_value: int = 0
    category: str = "autosomal"

    #: Kit-style allele number for the GRCh38 (or configured build) REF
    #: haplotype in this panel interval. Enables ``Δbp/step`` → absolute alleles for
    #: compound markers (period -1); curate per lab from STRBase / manufacturer tables.
    reference_ce: float | None = None
    #: Effective bp increment per +1 forensic repeat unit used with ``Δ bp`` versus
    #: the REF anchor length (defaults to tetranucleotide-style mapping).
    allele_bp_step: int = Field(default=4, ge=1)

    min_mapq: int | None = None
    min_mean_qual: int | None = None
    max_tr_len: int = 1000

    allow_triallelic: bool = False
    tri_balanced_thr: float | None = None
    #: Absolute read-count floor a 3rd candidate must clear (in addition to the
    #: fractional ``calling_thresh``) before it can promote a locus to triallelic
    #: or raise a ``mixture_suspected`` flag. Guards against ONT basecaller
    #: phantoms (2–4 reads ≈ 5–10% AF at 20–50× coverage). See known-bug #6.
    min_reads_third: int = Field(default=5, ge=0)
    ont_len_tolerance: int = 0
    stutter_overrides: dict[str, float] = Field(default_factory=dict)

    # "str" for normal TR markers; "amel" for Amelogenin sex typing.
    marker_type: str = "str"
    # Strand of the canonical ISFG motif relative to the reference (+/-).
    # Set to "-" for markers whose published motif is on the reverse strand
    # (e.g. D5S818, vWA, CSF1PO). The consensus is RC'd before ISFG compression.
    strand: str = "+"
    # AMELY region on chrY (AMEL markers only).
    y_chromosome: str | None = None
    y_ref_start: int | None = None
    y_ref_end: int | None = None

    @field_validator("strand")
    @classmethod
    def _strand_is_valid(cls, v: str) -> str:
        if v not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {v!r}")
        return v

    @field_validator("marker_type")
    @classmethod
    def _marker_type_is_known(cls, v: str) -> str:
        allowed = {"str", "amel"}
        if v not in allowed:
            raise ValueError(f"unknown marker_type {v!r}; expected one of {sorted(allowed)}")
        return v

    @field_validator("motif")
    @classmethod
    def _motif_is_nucleotide(cls, v: str) -> str:
        for m in v.split(","):
            if not m or any(c not in "ACGT" for c in m):
                raise ValueError(f"invalid motif {m!r}: must be non-empty [ACGT]+")
        return v

    @field_validator("category")
    @classmethod
    def _category_is_known(cls, v: str) -> str:
        allowed = {"autosomal", "y_chromosomal", "x_chromosomal", "mitochondrial"}
        if v not in allowed:
            raise ValueError(f"unknown category {v!r}; expected one of {sorted(allowed)}")
        return v

    def span(self) -> int:
        return self.ref_end - self.ref_start + 1


class Panel(BaseModel):
    """An ordered collection of systems with a version and provenance."""

    name: str
    version: str
    reference_build: str = "GRCh38"
    description: str | None = None
    systems: list[System]

    @field_validator("systems")
    @classmethod
    def _no_duplicate_names(cls, v: list[System]) -> list[System]:
        seen: set[str] = set()
        for s in v:
            if s.name in seen:
                raise ValueError(f"duplicate marker name in panel: {s.name!r}")
            seen.add(s.name)
        return v

    def by_name(self, name: str) -> System | None:
        for s in self.systems:
            if s.name == name:
                return s
        return None

    def by_chromosome(self) -> dict[str, list[System]]:
        out: dict[str, list[System]] = {}
        for s in self.systems:
            out.setdefault(s.chromosome, []).append(s)
        for v in out.values():
            v.sort(key=lambda x: x.ref_start)
        return out
