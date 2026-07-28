"""STRNaming-backed allele naming — the canonical source of the CE number.

The ISFG DNA Commission's current recommendations (Gettings et al. 2024,
FSI:Genetics 68:102946) name STRNaming as *the program* that must produce
bracketed repeat formatting, rather than describing a rule set to reimplement.
FRONTStr's hand-rolled :func:`~frontstr.interp.isfg.ce_from_brackets` counts raw
bracket units, and that count is provably wrong on six of the CODIS markers:
naming the GRCh38 reference window itself and comparing against the
authoritative ``ref_ce`` of STRNaming's range table shows vWA, D21S11, D2S1338,
D1S1656, D2S441 and DXS7132 all off. For the five compound (``period: -1``)
markers the offset is *allele-structure dependent*, so no per-marker
``corr_value`` can fix it — vWA and D21S11 were simply the two where it was
noticed first.

This module therefore delegates the number to STRNaming. It stays hermetic:
reference sequence comes from the committed slice cache built by
:mod:`frontstr.panel.seed_strnaming`, never from the network.

Why the reporting range is located by anchors, not coordinates
--------------------------------------------------------------

STRNaming's ranges hug the repeat array — TPOX's starts on the array's very
first base. An aligner is free to place a long allele's extra repeat units on
either side of that boundary, so slicing a read at the range's reference
coordinates silently drops them. Measured on HG00113: both TPOX alleles came
back as the 37 bp reference (CE8/CE8) when the truth is 9/11, and a 10-repeat
D8S1179 allele was named CE12.

So we locate the range by aligning the reference sequence *immediately flanking*
it against the cluster consensus and taking everything between the two hits.
That is indel-placement independent, and it is applied to the POA consensus
rather than to raw reads so ONT error does not move the anchors.
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache
from pathlib import Path

from frontstr.log import get_logger

log = get_logger(__name__)

#: Committed GRCh38 slice cache; see :mod:`frontstr.panel.seed_strnaming`.
CACHE_PATH = Path(__file__).parent / "data" / "strnaming_ranges.tsv"

#: bp of reference flank used to locate each end of a reporting range.
ANCHOR_BP = 30

#: Maximum edit distance, as a fraction of anchor length, for an anchor hit to
#: be trusted. edlib's infix mode always returns *some* location, so without
#: this a marker whose window does not actually contain the range would yield a
#: confident-looking nonsense slice.
MAX_ANCHOR_EDIT_FRACTION = 0.25

_CE_RE = re.compile(r"^CE(-?\d+(?:\.\d+)?)_")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

#: STRNaming's ``limit=`` option guards against sequences the Universal Analysis
#: Software truncates at a fixed length, returning a placeholder name when the
#: input is exactly that long. FRONTStr never truncates — it extracts a range
#: from a full-length consensus by anchors — so the guard only ever fires as a
#: false positive here. Left in place it makes DXS8378 unnameable even on the
#: reference sequence (``CE?_TODO_UAS_INCOMPLETE_SEQUENCE`` instead of ``CE10``).
_DROP_OPTIONS = ("limit",)

_RANGES_TSV_HEADER = "name\tchromosome\tstart\tend\tref_ce\toptions"


def reverse_complement(seq: str) -> str:
    """Reverse-complement a nucleotide sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


class NameStatus(StrEnum):
    """Why :meth:`StrNamer.name` did or did not produce a name."""

    OK = "ok"
    #: STRNaming defines no reported range for this marker (e.g. DYS393, AMEL).
    NO_RANGE = "no_range"
    #: One of the flanking anchors could not be placed in the consensus — the
    #: panel window is too narrow, or the consensus is too corrupt.
    ANCHOR_NOT_FOUND = "anchor_not_found"
    #: The consensus was empty (deletion / no data).
    EMPTY = "empty"
    #: STRNaming ran but emitted no ``CE<n>_`` prefix (no STR in range).
    NO_CE = "no_ce"
    #: STRNaming raised. The legacy CE path is used instead.
    ERROR = "error"


@dataclass(frozen=True, slots=True)
class NameResult:
    """Outcome of naming one allele consensus."""

    ce: float | None
    name: str
    status: NameStatus

    @property
    def ok(self) -> bool:
        return self.status is NameStatus.OK and self.ce is not None


@dataclass(frozen=True, slots=True)
class _Range:
    """One cached reporting range plus its extraction anchors."""

    name: str
    reverse: bool
    ref_ce: str
    left_anchor: str
    right_anchor: str


def _strip_dropped_options(options: str) -> str:
    """Remove options that do not apply to FRONTStr input (see ``_DROP_OPTIONS``)."""
    kept = [
        part
        for part in options.split(",")
        if part and not any(part.startswith(f"{name}=") for name in _DROP_OPTIONS)
    ]
    return ",".join(kept)


class StrNamer:
    """Names allele consensuses via STRNaming, fully offline.

    Construct with :meth:`from_cache`; the bundled instance is cached process-wide
    by :func:`default_namer` because building the reference structures costs
    real time and the result is immutable.
    """

    def __init__(self, store: object, ranges: dict[str, _Range]) -> None:
        self._store = store
        self._ranges = ranges

    @classmethod
    def from_cache(cls, path: Path | None = None) -> StrNamer:
        """Build a namer from the committed slice cache.

        Raises:
            InterpretationError: if strnaming is not installed or the cache is
                unreadable.
        """
        from frontstr.errors import InterpretationError

        path = path or CACHE_PATH
        try:
            from strnaming import classes
        except ImportError as exc:  # pragma: no cover - declared dependency
            raise InterpretationError(
                "strnaming is required for allele naming (pip install 'strnaming>=1.2,<1.3')"
            ) from exc

        try:
            text = path.read_text(encoding="utf-8")
        except OSError as exc:
            raise InterpretationError(f"Cannot read STRNaming range cache {path}: {exc}") from exc

        refseq_store = classes.ReferenceSequenceStore(autoload=False)
        tsv_rows: list[str] = []
        ranges: dict[str, _Range] = {}

        for line in text.splitlines()[1:]:
            if not line.strip():
                continue
            name, chromosome, start, end, ref_ce, options, slice_start, slice_seq = line.split("\t")
            start_i, end_i, slice_start_i = int(start), int(end), int(slice_start)
            refseq_store.add_refseq(chromosome, slice_start_i, slice_seq)

            options = _strip_dropped_options(options)
            tsv_rows.append("\t".join((name, chromosome, start, end, ref_ce, options)))

            lo, hi = min(start_i, end_i), max(start_i, end_i)
            # Offsets into slice_seq, whose first base is slice_start_i.
            left_end = lo - slice_start_i
            right_start = hi - slice_start_i + 1
            ranges[name] = _Range(
                name=name,
                reverse=end_i < start_i,
                ref_ce=ref_ce,
                left_anchor=slice_seq[max(0, left_end - ANCHOR_BP) : left_end],
                right_anchor=slice_seq[right_start : right_start + ANCHOR_BP],
            )

        structure_store = classes.ReferenceStructureStore(refseq_store=refseq_store)
        store = classes.ReportedRangeStore(structure_store=structure_store)
        store.load_from_tsv(
            io.StringIO(_RANGES_TSV_HEADER + "\n" + "\n".join(tsv_rows) + "\n"),
            load_structures=True,
        )
        return cls(store, ranges)

    def has_range(self, marker: str) -> bool:
        """True when STRNaming defines a reported range for ``marker``."""
        return marker in self._ranges

    def name(self, marker: str, consensus: str) -> NameResult:
        """Name one cluster consensus.

        Args:
            marker: Panel marker name.
            consensus: Cluster consensus over the **panel window**, in reference
                (forward) orientation — exactly what
                :class:`frontstr.evidence.cluster.Cluster` carries.

        Returns:
            A :class:`NameResult`. Any non-``OK`` status means the caller should
            fall back to the legacy CE path; this method never raises.
        """
        rng = self._ranges.get(marker)
        if rng is None:
            return NameResult(None, "", NameStatus.NO_RANGE)
        if not consensus:
            return NameResult(None, "", NameStatus.EMPTY)

        extracted = _extract_range(consensus, rng.left_anchor, rng.right_anchor)
        if extracted is None:
            log.debug("strnaming_anchor_failed", marker=marker, consensus_bp=len(consensus))
            return NameResult(None, "", NameStatus.ANCHOR_NOT_FOUND)

        query = reverse_complement(extracted) if rng.reverse else extracted
        try:
            name = self._store.get_range(marker).get_name(query)  # type: ignore[attr-defined]
        except Exception as exc:
            log.warning("strnaming_failed", marker=marker, error=str(exc))
            return NameResult(None, "", NameStatus.ERROR)

        match = _CE_RE.match(name)
        if match is None:
            return NameResult(None, name, NameStatus.NO_CE)
        return NameResult(float(match.group(1)), name, NameStatus.OK)


@lru_cache(maxsize=1)
def default_namer() -> StrNamer | None:
    """The process-wide namer over the bundled cache, or ``None`` if unavailable.

    Returning ``None`` rather than raising keeps naming an enhancement: a broken
    or missing cache degrades the run to the legacy CE path instead of killing
    it.
    """
    try:
        return StrNamer.from_cache()
    except Exception as exc:
        log.warning("strnaming_unavailable", error=str(exc))
        return None


def _locate(anchor: str, target: str) -> tuple[int, int] | None:
    """Best infix location of ``anchor`` within ``target`` as ``(start, end)``."""
    import edlib

    if not anchor or not target:
        return None
    result = edlib.align(anchor, target, mode="HW", task="locations")
    locations = result.get("locations") or []
    if not locations:
        return None
    start, end = locations[0]
    if start is None or end is None:
        return None
    return start, end + 1


def _extract_range(consensus: str, left_anchor: str, right_anchor: str) -> str | None:
    """Slice the reporting range out of ``consensus`` using its flanking anchors.

    Returns ``None`` when either anchor cannot be placed confidently, or when
    the hits are ordered such that no range lies between them.
    """
    import edlib

    left = _locate(left_anchor, consensus)
    right = _locate(right_anchor, consensus)
    if left is None or right is None:
        return None

    for anchor, hit in ((left_anchor, left), (right_anchor, right)):
        distance = edlib.align(anchor, consensus[hit[0] : hit[1]], mode="NW", task="distance")
        if distance["editDistance"] > MAX_ANCHOR_EDIT_FRACTION * len(anchor):
            return None

    start, end = left[1], right[0]
    if end <= start:
        return None
    return consensus[start:end]
