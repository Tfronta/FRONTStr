"""Per-locus read pileup: extracts the TR subsequence from each spanning read.

This is **Layer 2** of FRONTStr: integer per-allele coverage straight from the
BAM, independent of the caller. See plan-longtr-improved.md §5.2.

Coordinate convention
---------------------

Internally we use **0-based half-open** intervals (``[start, end)``), which is
what pysam and BED conventionally use. Callers that come from forensic panel
YAML (1-based inclusive) must convert:

    start_0based = system.ref_start - 1
    end_0based   = system.ref_end          # already exclusive of the next base

Strand handling
---------------

Per the SAM v1 specification, the ``SEQ`` field of a reverse-strand read is
already reverse-complemented to forward orientation. ``pysam.query_sequence``
returns the SAM ``SEQ`` field verbatim, so ``query_sequence[a:b]`` is always
in **reference (forward)** orientation regardless of read strand. We record
strand in the :class:`Observation` for QC but the sequence is canonical.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any

from frontstr.errors import EvidenceError


@dataclass(slots=True, frozen=True)
class Observation:
    """A single read's contribution at one locus.

    Attributes:
        read_id: BAM ``QNAME`` (used for auditing / drill-down).
        sequence: TR subsequence in reference orientation, uppercase.
        hp: Haplotype tag (``1``, ``2``) or ``None`` if absent.
        mean_qual: Mean Phred quality of the TR subsequence.
        strand: ``"+"`` if mapped to forward strand, ``"-"`` for reverse.
        flank_left_ok: True if the left flank anchor was clean (no clip/indel).
        flank_right_ok: Same for the right flank anchor.
    """

    read_id: str
    sequence: str
    hp: int | None
    mean_qual: float
    strand: str
    flank_left_ok: bool
    flank_right_ok: bool


_DEFAULT_FLANK_ANCHOR = 20
_DEFAULT_FETCH_MARGIN = 200


class RejectReason(StrEnum):
    """Why a fetched read did not become an :class:`Observation`.

    These used to be invisible: :func:`_read_to_observation` returned ``None``
    and the caller saw only the survivors, so "35 reads at this locus" could not
    be distinguished from "35 of 41, and here is what happened to the other 6".
    For a caller whose output is used forensically that distinction is the
    difference between a coverage number and an auditable one.
    """

    NOT_PRIMARY = "not a primary alignment"
    LOW_MAPQ = "MAPQ below threshold"
    NO_ALIGNMENT_END = "no alignment end (unmapped or malformed)"
    LEFT_FLANK_SHORT = "does not reach the left flank anchor"
    RIGHT_FLANK_SHORT = "does not reach the right flank anchor"
    NO_SEQUENCE = "read carries no sequence"
    WINDOW_UNLOCATABLE = "window could not be located in the read"


@dataclass(slots=True)
class PileupCounts:
    """Read accounting for one locus: how many came in, what left, and why.

    Passed into :func:`pileup_locus` to be filled. Optional so the hot path
    stays allocation-free when nobody is tracing.
    """

    fetched: int = 0
    kept: int = 0
    rejected: Counter[RejectReason] = field(default_factory=Counter)

    @property
    def n_rejected(self) -> int:
        return sum(self.rejected.values())

    def reasons(self) -> list[tuple[RejectReason, int]]:
        """Rejection reasons, most common first, for display."""
        return sorted(self.rejected.items(), key=lambda kv: (-kv[1], kv[0].value))


def pileup_locus(
    bam_path: Path,
    chrom: str,
    start: int,
    end: int,
    *,
    min_mapq: int = 20,
    flank_anchor: int = _DEFAULT_FLANK_ANCHOR,
    fetch_margin: int = _DEFAULT_FETCH_MARGIN,
    reference_fasta: Path | None = None,
    counts: PileupCounts | None = None,
) -> list[Observation]:
    """Extract one :class:`Observation` per read that fully spans ``[start, end)``.

    Args:
        bam_path: Indexed BAM (``.bai`` or ``.csi`` next to it), or an indexed
            CRAM (``.crai`` next to it).
        chrom: Chromosome name as it appears in the BAM ``@SQ`` header.
        start: 0-based inclusive start of the TR.
        end: 0-based exclusive end of the TR.
        min_mapq: Drop reads with MAPQ below this threshold.
        flank_anchor: Required bp of clean reference-aligned flank on each side
            (drops reads whose alignment edge falls inside this window).
        fetch_margin: How far around the locus to ``fetch`` reads. Reads whose
            alignment extends beyond ``start - fetch_margin`` to ``end +
            fetch_margin`` are inspected; the actual gating is by alignment
            boundary checks below.
        reference_fasta: Path to the reference FASTA (with ``.fai`` index).
            Required for CRAM input; ignored for BAM.
        counts: Optional :class:`PileupCounts` to fill with the read funnel —
            how many were fetched, how many survived and why the rest did not.
            Left ``None`` the accounting is not computed at all.

    Returns:
        A list of :class:`Observation`, possibly empty. Order matches the BAM
        iteration order (which is coordinate-sorted for sorted BAMs).

    Raises:
        EvidenceError: If the BAM/CRAM cannot be opened or fetch fails.
    """
    if start < 0 or end <= start:
        raise EvidenceError(f"Invalid locus window: start={start}, end={end}")
    if not bam_path.exists():
        raise EvidenceError(f"BAM not found: {bam_path}")

    is_cram = bam_path.suffix.lower() == ".cram"
    if is_cram and reference_fasta is None:
        raise EvidenceError(
            "CRAM file requires a reference FASTA — "
            "pass reference_fasta= to pileup_locus (or --reference on the CLI)"
        )

    try:
        import pysam
    except ImportError as exc:
        raise EvidenceError("pysam is required for pileup_locus") from exc

    open_kwargs: dict[str, Any] = {"reference_filename": str(reference_fasta)} if is_cram else {}
    try:
        bam = pysam.AlignmentFile(str(bam_path), "rc" if is_cram else "rb", **open_kwargs)
    except (ValueError, OSError) as exc:
        raise EvidenceError(f"Cannot open BAM {bam_path}: {exc}") from exc

    obs: list[Observation] = []
    fetch_start = max(0, start - fetch_margin)
    fetch_end = end + fetch_margin
    try:
        for read in bam.fetch(chrom, fetch_start, fetch_end):
            if counts is not None:
                counts.fetched += 1
            o, reason = _read_to_observation(read, start, end, min_mapq, flank_anchor)
            if o is not None:
                obs.append(o)
                if counts is not None:
                    counts.kept += 1
            elif counts is not None and reason is not None:
                counts.rejected[reason] += 1
    except (ValueError, OSError) as exc:
        raise EvidenceError(f"Fetch failed on {chrom}:{start}-{end}: {exc}") from exc
    finally:
        bam.close()
    return obs


def _read_to_observation(
    read: object,  # pysam.AlignedSegment — typed as object to keep pysam optional
    start: int,
    end: int,
    min_mapq: int,
    flank_anchor: int,
) -> tuple[Observation | None, RejectReason | None]:
    """Decide whether ``read`` contributes an :class:`Observation` and build it.

    Returns ``(observation, None)`` on success and ``(None, reason)`` otherwise.
    The reason is returned rather than swallowed so a locus can account for
    every read it fetched — see :class:`PileupCounts`.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:  # type: ignore[attr-defined]
        return None, RejectReason.NOT_PRIMARY
    if read.mapping_quality < min_mapq:  # type: ignore[attr-defined]
        return None, RejectReason.LOW_MAPQ
    ref_start: int = read.reference_start  # type: ignore[attr-defined]
    ref_end: int | None = read.reference_end  # type: ignore[attr-defined]
    if ref_end is None:
        return None, RejectReason.NO_ALIGNMENT_END
    if ref_start > start - flank_anchor:
        return None, RejectReason.LEFT_FLANK_SHORT
    if ref_end < end + flank_anchor:
        return None, RejectReason.RIGHT_FLANK_SHORT

    qseq: str | None = read.query_sequence  # type: ignore[attr-defined]
    if not qseq:
        return None, RejectReason.NO_SEQUENCE

    qstart, qend = _locate_window(read, start, end)
    if qstart is None or qend is None or qend <= qstart:
        return None, RejectReason.WINDOW_UNLOCATABLE

    sequence = qseq[qstart:qend].upper()

    hp: int | None = None
    try:
        if read.has_tag("HP"):  # type: ignore[attr-defined]
            hp_val = read.get_tag("HP")  # type: ignore[attr-defined]
            hp = int(hp_val) if hp_val in (1, 2, "1", "2") else None
    except (ValueError, KeyError):
        hp = None

    mean_qual = _mean_quality(read, qstart, qend)
    strand = "-" if read.is_reverse else "+"  # type: ignore[attr-defined]

    return (
        Observation(
            read_id=read.query_name or "?",  # type: ignore[attr-defined]
            sequence=sequence,
            hp=hp,
            mean_qual=mean_qual,
            strand=strand,
            flank_left_ok=True,
            flank_right_ok=True,
        ),
        None,
    )


def _locate_window(
    read: object,  # pysam.AlignedSegment
    start: int,
    end: int,
) -> tuple[int | None, int | None]:
    """Translate the reference window ``[start, end)`` into query coordinates.

    Returns ``(qstart, qend)`` such that ``query_sequence[qstart:qend]`` is the
    portion of the read aligned to the requested reference window. Either may
    be ``None`` if the alignment does not cover the requested position.

    Handles insertions (extra read bases inside the window) and deletions
    (missing reference positions inside the window) naturally because
    ``get_aligned_pairs`` reports indel positions with ``None`` on one side.

    Deletion-at-boundary handling:
    - If a deletion covers ``start``, qstart is set to the first non-None qpos
      at any rpos >= start (i.e., the first read base after the deleted region).
    - If a deletion covers ``end``, qend is set to last_qpos + 1 (exclusive
      slice index after the last read base that was inside the window).
    """
    qstart: int | None = None
    qend: int | None = None
    last_qpos: int | None = None  # tracks last non-None qpos for deletion-at-end
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):  # type: ignore[attr-defined]
        if rpos is None:
            continue
        if qpos is not None:
            last_qpos = qpos
        # Check end boundary first so we break as soon as the window is passed.
        if rpos >= end:
            if qpos is not None and rpos == end:
                qend = qpos
            else:
                # Deletion at or spanning end: use position after last matched base.
                qend = last_qpos + 1 if last_qpos is not None else None
            break
        # Accept the first non-None qpos at any rpos >= start so that a deletion
        # covering the exact start position does not silently lose the read.
        if qstart is None and rpos >= start and qpos is not None:
            qstart = qpos
    return qstart, qend


def _mean_quality(read: object, qstart: int, qend: int) -> float:
    """Mean Phred quality across ``query_qualities[qstart:qend]``."""
    quals = read.query_qualities  # type: ignore[attr-defined]
    if quals is None:
        return 0.0
    window = list(quals)[qstart:qend]
    if not window:
        return 0.0
    return float(sum(window)) / len(window)
