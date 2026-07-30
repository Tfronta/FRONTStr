"""Allele catalog: sequence to ISFG/CE lookup table.

This is the file-based MVP analogue of ``allele_catalog`` from
plan-longtr-improved.md §3.1. Instead of a Postgres table we ship a committable
JSON document (``examples/catalogs/strseq_*.json``) loaded into an in-memory
:class:`AlleleCatalog`.

A :class:`CatalogEntry` is one *core* allele sequence (flanks stripped) with its
forensically curated ISFG string, CE number, and ISFG iso-allele suffix
(``a``/``b``/``c`` …). The catalog lets the interpretation layer report the
**published** ISFG/suffix for a known sequence rather than a freshly computed one,
and to distinguish iso-alleles that share a CE but differ in internal structure
(e.g. D3S1358 ``14a`` vs ``14b``).

Lookup is keyed by marker name; within a marker, candidates are pre-filtered by
length and refined by edit distance in :mod:`frontstr.interp.catalog`.
"""

from __future__ import annotations

import json
from pathlib import Path

from pydantic import BaseModel, Field, field_validator

from frontstr.errors import PanelError


class CatalogEntry(BaseModel):
    """One curated core-sequence → ISFG/CE record for a single marker."""

    marker: str
    #: Core repeat sequence with flanks stripped, uppercase [ACGT].
    sequence: str
    #: Curated ISFG bracketed string, e.g. ``"[TCTG]2 [TCTA]11"``.
    isfg: str
    #: Forensic CE allele number (e.g. ``14``, ``9.3``); ``None`` if undefined.
    ce: float | None = None
    #: ISFG iso-allele suffix distinguishing same-CE structures ("a", "b", …).
    suffix: str | None = None
    #: Provenance: ``"STRSeq"``, ``"lab-internal"``, ``"manual"``.
    source: str = "manual"
    #: GenBank accession when ``source == "STRSeq"``.
    accession: str | None = None
    #: Optional flanks for realign cross-checks (50 bp each side).
    flank_left: str | None = None
    flank_right: str | None = None

    @field_validator("sequence")
    @classmethod
    def _sequence_is_nucleotide(cls, v: str) -> str:
        up = v.upper()
        if not up or any(c not in "ACGT" for c in up):
            raise ValueError(f"catalog sequence must be non-empty [ACGT]+, got {v!r}")
        return up


class AlleleCatalog(BaseModel):
    """A versioned collection of :class:`CatalogEntry`, indexed by marker."""

    version: str
    reference_build: str = "GRCh38"
    source: str = "STRSeq"
    description: str | None = None
    entries: list[CatalogEntry] = Field(default_factory=list)

    def by_marker(self, marker: str) -> list[CatalogEntry]:
        """All catalog entries for ``marker`` (empty list if none)."""
        return self._index.get(marker, [])

    @property
    def _index(self) -> dict[str, list[CatalogEntry]]:
        # Built lazily and cached on the instance (pydantic private-attr free).
        idx = self.__dict__.get("_marker_index")
        if idx is None:
            idx = {}
            for e in self.entries:
                idx.setdefault(e.marker, []).append(e)
            self.__dict__["_marker_index"] = idx
        return idx

    def markers(self) -> list[str]:
        return sorted(self._index)


def load_catalog(path: Path) -> AlleleCatalog:
    """Load and validate an :class:`AlleleCatalog` from a JSON file."""
    if not path.exists():
        raise PanelError(f"catalog file not found: {path}")
    try:
        raw = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise PanelError(f"catalog JSON parse error in {path}: {exc}") from exc
    try:
        return AlleleCatalog.model_validate(raw)
    except Exception as exc:  # pydantic ValidationError
        raise PanelError(f"invalid catalog in {path}: {exc}") from exc


def dump_catalog(catalog: AlleleCatalog, path: Path, *, indent: int = 2) -> None:
    """Serialize ``catalog`` to JSON at ``path`` (deterministic key order)."""
    data = catalog.model_dump(exclude_none=True)
    path.write_text(json.dumps(data, indent=indent, sort_keys=False) + "\n")
