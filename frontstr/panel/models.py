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
    ont_len_tolerance: int = 0
    stutter_overrides: dict[str, float] = Field(default_factory=dict)

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
