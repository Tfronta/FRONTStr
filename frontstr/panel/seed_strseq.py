"""Build an :class:`AlleleCatalog` from STRSeq records (Phase 2.1 importer).

STRSeq is the NIST/ISFG curated STR sequence resource published on GenBank as
two BioProjects (same source STRspy 2.0 seeds from):

- ``PRJNA380345`` — autosomal STR alleles
- ``PRJNA380347`` — Y-STR alleles

This module separates the two concerns so the assembly logic is unit-testable
without network access:

1. :func:`build_catalog` — pure function: a list of parsed records (marker,
   core sequence, CE, suffix, accession) → validated :class:`AlleleCatalog`.
   This is what the committable ``examples/catalogs/strseq_*.json`` is produced
   with, and it is fully covered by tests.
2. :func:`fetch_strseq_records` — the NCBI Entrez fetch + GenBank feature
   parse. This is the network/curation step; it is intentionally a thin,
   clearly-marked boundary because the flank-trimming rules are marker-specific
   and must be curated against each panel's windows before the output is
   forensically usable.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from frontstr.panel.catalog import AlleleCatalog, CatalogEntry, dump_catalog

#: STRSeq GenBank BioProjects (see module docstring).
STRSEQ_BIOPROJECTS = {
    "autosomal": "PRJNA380345",
    "y_chromosomal": "PRJNA380347",
}


@dataclass(slots=True)
class StrSeqRecord:
    """One parsed STRSeq allele record (flanks already trimmed to the core)."""

    marker: str
    sequence: str
    isfg: str
    ce: float | None = None
    suffix: str | None = None
    accession: str | None = None
    flank_left: str | None = None
    flank_right: str | None = None


def build_catalog(
    records: Iterable[StrSeqRecord],
    *,
    version: str,
    reference_build: str = "GRCh38",
    source: str = "STRSeq",
    description: str | None = None,
) -> AlleleCatalog:
    """Assemble a validated :class:`AlleleCatalog` from parsed STRSeq records.

    Pure and deterministic — the unit-tested half of the importer. Duplicate
    ``(marker, sequence)`` pairs are collapsed, keeping the first occurrence.
    """
    seen: set[tuple[str, str]] = set()
    entries: list[CatalogEntry] = []
    for r in records:
        key = (r.marker, r.sequence.upper())
        if key in seen:
            continue
        seen.add(key)
        entries.append(
            CatalogEntry(
                marker=r.marker,
                sequence=r.sequence,
                isfg=r.isfg,
                ce=r.ce,
                suffix=r.suffix,
                source=source,
                accession=r.accession,
                flank_left=r.flank_left,
                flank_right=r.flank_right,
            )
        )
    return AlleleCatalog(
        version=version,
        reference_build=reference_build,
        source=source,
        description=description,
        entries=entries,
    )


def write_catalog(catalog: AlleleCatalog, out_path: Path) -> Path:
    """Dump ``catalog`` to ``out_path`` as JSON and return the path."""
    dump_catalog(catalog, out_path)
    return out_path


def fetch_strseq_records(
    bioproject: str,
    *,
    email: str,
) -> list[StrSeqRecord]:  # pragma: no cover - network/curation boundary
    """Fetch + parse STRSeq GenBank records for a BioProject.

    NOT YET IMPLEMENTED. This is the network/curation step of the importer.
    To complete it:

    1. Use Bio.Entrez (esearch over ``bioproject`` → efetch GenBank flat files)
       with the caller-supplied ``email`` (NCBI requires it).
    2. For each record, read the STR repeat feature, derive the core sequence,
       and **trim flanks to match the panel extraction window** (this is the
       marker-specific curation that must be reviewed before use).
    3. Compute or carry the ISFG string + CE + iso-allele suffix.

    Return the parsed records and pass them to :func:`build_catalog`. Keep this
    boundary thin so the assembly logic above stays offline-testable.
    """
    raise NotImplementedError(
        "STRSeq GenBank fetch is the network/curation step — see docstring. "
        f"BioProject {bioproject!r} (known: {sorted(STRSEQ_BIOPROJECTS.values())})."
    )
