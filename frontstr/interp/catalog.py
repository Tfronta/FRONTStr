"""Annotate alleles with a curated allele catalog (Phase 2.1 / plan §6.1).

``annotate_with_catalog`` takes one classified :class:`Allele`, looks up the
marker's :class:`CatalogEntry` set, and assigns the **published** ISFG / CE /
iso-allele suffix when the cluster consensus matches a known sequence.

Decision (mirrors plan-longtr-improved.md §6.1):

- edit distance ``0`` → exact match: adopt the catalog ISFG + suffix
  (and CE when the live value is missing).
- ``0 < distance <= eps`` → approximate match: adopt the catalog ISFG, suffix
  gets a trailing ``*`` to flag the imperfect hit.
- otherwise → no catalog hit: leave the live-computed ISFG/CE untouched.

The function mutates ``allele`` in place and returns a :class:`CatalogMatch`
describing what happened (useful for tests, reports, and audit).
"""

from __future__ import annotations

from dataclasses import dataclass

from frontstr.interp.isfg import _rc
from frontstr.interp.models import Allele, IsoAllele
from frontstr.interp.stutter import MotifRun, find_motif_runs
from frontstr.panel.catalog import AlleleCatalog, CatalogEntry
from frontstr.panel.models import System

#: Max edit distance still accepted as an (approximate) catalog hit.
DEFAULT_EXACT_EPS = 2
#: Candidates whose length differs from the consensus by more than this are
#: skipped before the (costlier) edit-distance refinement.
DEFAULT_LEN_WINDOW = 6
#: Non-motif chars between two motif runs above which they belong to different
#: blocks. Matches ``calibrate._GAP_THRESHOLD``: intra-allele spacers (ATG, ta,
#: tccata …) are ≤6 chars, so 12 cleanly separates the repeat core from chance
#: motif hits in the ±100bp flanks of the extraction window.
DEFAULT_GAP_THRESHOLD = 12


def extract_repeat_core(
    sequence: str,
    motifs: list[str],
    *,
    strand: str = "+",
    gap_threshold: int = DEFAULT_GAP_THRESHOLD,
) -> str:
    """Reduce a flank-inclusive consensus to its repeat-core substring.

    FRONTStr clusters carry the full ±100bp extraction window; catalog entries
    store only the repeat core. This finds the motif runs, groups runs separated
    by ≤ ``gap_threshold`` non-motif chars, and returns the nucleotide substring
    spanning the unit-richest group (internal spacers preserved, flanks and
    chance flank motif hits dropped). Returns the input unchanged if no run is
    found. Reverse-strand markers are canonicalised first.
    """
    if not sequence or not motifs:
        return sequence
    seq = _rc(sequence) if strand == "-" else sequence
    runs = find_motif_runs(seq, motifs)
    if not runs:
        return seq

    groups: list[list[MotifRun]] = [[runs[0]]]
    for run in runs[1:]:
        if run.start - groups[-1][-1].end <= gap_threshold:
            groups[-1].append(run)
        else:
            groups.append([run])

    core = max(groups, key=lambda g: sum(r.length_bp for r in g))
    return seq[core[0].start : core[-1].end]


@dataclass(slots=True)
class CatalogMatch:
    """Outcome of a catalog lookup for one allele."""

    matched: bool
    exact: bool
    distance: int | None
    entry: CatalogEntry | None


def _edit_distance(a: str, b: str) -> int:
    import edlib

    return int(edlib.align(a, b, task="distance", mode="NW")["editDistance"])


def best_catalog_hit(
    sequence: str,
    candidates: list[CatalogEntry],
    *,
    len_window: int = DEFAULT_LEN_WINDOW,
) -> tuple[CatalogEntry | None, int | None]:
    """Return the lowest-edit-distance catalog entry for ``sequence``.

    Candidates farther than ``len_window`` bp in length are skipped. Returns
    ``(None, None)`` when no candidate survives the length filter.
    """
    if not sequence or not candidates:
        return None, None
    best: CatalogEntry | None = None
    best_d: int | None = None
    for c in candidates:
        if abs(len(c.sequence) - len(sequence)) > len_window:
            continue
        d = _edit_distance(sequence, c.sequence)
        if best_d is None or d < best_d:
            best, best_d = c, d
    return best, best_d


def annotate_with_catalog(
    allele: Allele,
    system: System,
    catalog: AlleleCatalog,
    *,
    exact_eps: int = DEFAULT_EXACT_EPS,
    len_window: int = DEFAULT_LEN_WINDOW,
) -> CatalogMatch:
    """Annotate ``allele`` from ``catalog`` in place; return the match outcome."""
    if not allele.consensus or allele.is_deletion:
        return CatalogMatch(matched=False, exact=False, distance=None, entry=None)

    # Catalog sequences are repeat cores; the cluster consensus is flank-inclusive.
    motifs = [m for m in system.motif.split(",") if m]
    seq = extract_repeat_core(allele.consensus, motifs, strand=system.strand)

    candidates = catalog.by_marker(system.name)
    entry, dist = best_catalog_hit(seq, candidates, len_window=len_window)

    if entry is None or dist is None or dist > exact_eps:
        # No usable hit — keep the live-computed ISFG/CE.
        return CatalogMatch(matched=False, exact=False, distance=dist, entry=None)

    allele.isfg = entry.isfg
    if allele.ce is None and entry.ce is not None:
        allele.ce = entry.ce

    if dist == 0:
        allele.iso = IsoAllele(
            suffix=entry.suffix, match_type="exact", distance=0, source=entry.source
        )
        return CatalogMatch(matched=True, exact=True, distance=0, entry=entry)

    # Approximate hit: clean suffix, match_type records the inexactness.
    allele.iso = IsoAllele(
        suffix=entry.suffix, match_type="approx", distance=dist, source=entry.source
    )
    return CatalogMatch(matched=True, exact=False, distance=dist, entry=entry)


def annotate_alleles(
    alleles: list[Allele],
    system: System,
    catalog: AlleleCatalog | None,
    **kwargs: object,
) -> None:
    """Annotate every allele in ``alleles`` against ``catalog`` (no-op if None)."""
    if catalog is None:
        return
    for a in alleles:
        annotate_with_catalog(a, system, catalog, **kwargs)  # type: ignore[arg-type]
