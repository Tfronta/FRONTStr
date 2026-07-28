"""Build the offline GRCh38 slice cache that :mod:`frontstr.interp.naming` reads.

STRNaming needs reference sequence to work out each reported range's repeat
structure, and by default it fetches that from the Ensembl REST API into
``~/.strnaming-cache``. A forensic caller must not put a network round-trip in
the middle of allele naming, and ``frontstr interpret`` deliberately does not
require a reference FASTA for BAM input, so neither of those routes is
acceptable here.

Instead we commit the few kilobases that actually matter. For every marker in
the panel that STRNaming defines a reported range for, this module extracts the
range plus :data:`SLICE_PAD` bp on each side and writes it to a TSV that ships
inside the package. At run time :mod:`frontstr.interp.naming` seeds a
``ReferenceSequenceStore(autoload=False)`` from that file, so naming is
hermetic: no network, no FASTA, byte-identical results on every machine.

This is a **calibration-time** step, like :mod:`frontstr.panel.calibrate`. It
needs an indexed GRCh38 FASTA and is only re-run when the panel gains markers
or STRNaming ships new ranges.
"""

from __future__ import annotations

import pkgutil
from dataclasses import dataclass
from pathlib import Path

from frontstr.errors import PanelError
from frontstr.panel.models import Panel

#: bp of reference kept on each side of a reported range. STRNaming needs the
#: repeat structures that overlap the range, which can extend past it, plus
#: :data:`frontstr.interp.naming.ANCHOR_BP` for the extraction anchors.
#: Empirically 100 bp suffices for every CODIS marker; 200 is the safety margin
#: and still leaves the whole cache around 15 kB.
SLICE_PAD = 200

#: STRNaming's own per-marker range table (name, chromosome, start, end,
#: ref_ce, options). ``start > end`` encodes a minus-strand range.
STRNAMING_RANGES_RESOURCE = "data/ranges_uas-frr.txt"

CACHE_HEADER = (
    "name",
    "chromosome",
    "start",
    "end",
    "ref_ce",
    "options",
    "slice_start",
    "slice_seq",
)


@dataclass(frozen=True, slots=True)
class RangeRow:
    """One marker's reported range plus the reference slice backing it."""

    name: str
    chromosome: str
    start: int
    end: int
    ref_ce: str
    options: str
    #: 1-based genome position of the first base of :attr:`slice_seq`.
    slice_start: int
    slice_seq: str

    @property
    def lo(self) -> int:
        """Lowest 1-based coordinate of the range (inclusive)."""
        return min(self.start, self.end)

    @property
    def hi(self) -> int:
        """Highest 1-based coordinate of the range (inclusive)."""
        return max(self.start, self.end)


def load_strnaming_ranges() -> dict[str, tuple[str, int, int, str, str]]:
    """Parse STRNaming's bundled range table.

    Returns:
        ``{marker: (chromosome, start, end, ref_ce, options)}`` with start/end
        exactly as STRNaming writes them (start > end means minus strand).
    """
    raw = pkgutil.get_data("strnaming", STRNAMING_RANGES_RESOURCE)
    if raw is None:  # pragma: no cover - packaging accident
        raise PanelError(f"strnaming is missing {STRNAMING_RANGES_RESOURCE}")
    out: dict[str, tuple[str, int, int, str, str]] = {}
    for line in raw.decode().splitlines()[1:]:
        if not line.strip():
            continue
        f = line.split("\t")
        out[f[0]] = (f[1], int(f[2]), int(f[3]), f[4], f[5] if len(f) > 5 else "")
    return out


def build_rows(panel: Panel, fetch_seq: object) -> tuple[list[RangeRow], list[str]]:
    """Assemble cache rows for every panel marker STRNaming defines a range for.

    Args:
        panel: The panel whose markers should be cached.
        fetch_seq: Callable ``(chromosome, start, end) -> str`` returning
            uppercase reference sequence for a 1-based inclusive interval.
            Injected so the assembly logic stays testable without a FASTA.

    Returns:
        ``(rows, skipped)`` where ``skipped`` names the panel markers STRNaming
        has no reported range for (they keep the legacy CE path).
    """
    defined = load_strnaming_ranges()
    rows: list[RangeRow] = []
    skipped: list[str] = []
    for system in panel.systems:
        entry = defined.get(system.name)
        if entry is None:
            skipped.append(system.name)
            continue
        chromosome, start, end, ref_ce, options = entry
        lo, hi = min(start, end), max(start, end)
        slice_start = max(1, lo - SLICE_PAD)
        seq = fetch_seq(chromosome, slice_start, hi + SLICE_PAD)  # type: ignore[operator]
        rows.append(
            RangeRow(
                name=system.name,
                chromosome=chromosome,
                start=start,
                end=end,
                ref_ce=ref_ce,
                options=options,
                slice_start=slice_start,
                slice_seq=seq.upper(),
            )
        )
    return rows, skipped


def dump_cache(rows: list[RangeRow], path: Path) -> None:
    """Write cache rows as TSV (one line per marker)."""
    lines = ["\t".join(CACHE_HEADER)]
    lines.extend(
        "\t".join(
            (
                r.name,
                r.chromosome,
                str(r.start),
                str(r.end),
                r.ref_ce,
                r.options,
                str(r.slice_start),
                r.slice_seq,
            )
        )
        for r in sorted(rows, key=lambda r: r.name)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fasta_fetcher(reference_fasta: Path) -> object:
    """Return a ``(chromosome, start, end) -> str`` reader over an indexed FASTA.

    Accepts STRNaming's bare chromosome names (``12``, ``X``) and maps them onto
    whichever naming the FASTA uses (``chr12`` or ``12``).
    """
    import pysam

    fa = pysam.FastaFile(str(reference_fasta))
    available = set(fa.references)

    def fetch(chromosome: str, start: int, end: int) -> str:
        for candidate in (f"chr{chromosome}", chromosome):
            if candidate in available:
                return fa.fetch(candidate, start - 1, end).upper()
        raise PanelError(f"Reference {reference_fasta} has no contig for chromosome {chromosome!r}")

    return fetch
