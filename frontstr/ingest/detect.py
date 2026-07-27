"""Input format detection by magic bytes.

Forensic uploads frequently arrive with non-canonical filenames (e.g.
``S001_final_v3.dat``) so extension-based detection is unsafe. This module
inspects file content and answers three questions:

1. What is the format? (``fastq``, ``fasta``, ``bam``, ``cram``)
2. Is it gzipped (BGZF or plain gzip)?
3. If it is a BAM/CRAM, is it actually aligned, or just unaligned reads (uBAM)?

The function performs **bounded reads** (a few kB) so it is safe to call on
huge inputs without buffering them in memory.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from frontstr.errors import IngestError


class InputKind(StrEnum):
    """Recognised input format categories."""

    FASTQ = "fastq"
    FASTA = "fasta"
    BAM = "bam"
    UBAM = "ubam"  # unaligned BAM
    CRAM = "cram"


@dataclass(frozen=True, slots=True)
class InputInfo:
    """Outcome of :func:`detect_input`.

    Attributes:
        kind: Format category.
        gzipped: True if the on-disk file is bgzf/gzip-wrapped.
                 BAMs are always gzipped (bgzf) so this is True for them.
        aligned: For BAM/CRAM, whether at least one alignment record is mapped.
                 None for FASTA/FASTQ.
        size_bytes: Size of the file on disk.
    """

    kind: InputKind
    gzipped: bool
    aligned: bool | None
    size_bytes: int


_BAM_MAGIC = b"BAM\x01"
_CRAM_MAGIC = b"CRAM"
_GZIP_MAGIC = b"\x1f\x8b"


def detect_input(path: Path) -> InputInfo:
    """Detect the format of ``path`` by inspecting its magic bytes.

    Args:
        path: Path to an input file.

    Returns:
        :class:`InputInfo` describing the format.

    Raises:
        IngestError: If the file is empty, unreadable, or of an unsupported format.
    """
    if not path.exists():
        raise IngestError(f"Input does not exist: {path}")
    size = path.stat().st_size
    if size == 0:
        raise IngestError(f"Input is empty: {path}")

    try:
        with path.open("rb") as fh:
            head = fh.read(8)
    except OSError as exc:
        raise IngestError(f"Cannot read input {path}: {exc}") from exc

    if head[:4] == _CRAM_MAGIC:
        aligned = _bam_or_cram_has_alignment(path)
        return InputInfo(kind=InputKind.CRAM, gzipped=False, aligned=aligned, size_bytes=size)

    if head[:2] == _GZIP_MAGIC:
        return _classify_gzipped(path, size)

    if head[:1] == b"@":
        return InputInfo(kind=InputKind.FASTQ, gzipped=False, aligned=None, size_bytes=size)
    if head[:1] == b">":
        return InputInfo(kind=InputKind.FASTA, gzipped=False, aligned=None, size_bytes=size)

    raise IngestError(
        f"Unrecognised input format for {path}: head bytes {head[:8]!r} not FASTA/FASTQ/BAM/CRAM"
    )


def _classify_gzipped(path: Path, size: int) -> InputInfo:
    """Decide whether a gzipped file is BAM or fastq/fasta.gz."""
    try:
        with gzip.open(path, "rb") as fh:
            inner_head = fh.read(8)
    except OSError as exc:
        raise IngestError(f"Cannot read gzipped input {path}: {exc}") from exc

    if inner_head[:4] == _BAM_MAGIC:
        aligned = _bam_or_cram_has_alignment(path)
        kind = InputKind.BAM if aligned else InputKind.UBAM
        return InputInfo(kind=kind, gzipped=True, aligned=aligned, size_bytes=size)

    if inner_head[:1] == b"@":
        return InputInfo(kind=InputKind.FASTQ, gzipped=True, aligned=None, size_bytes=size)
    if inner_head[:1] == b">":
        return InputInfo(kind=InputKind.FASTA, gzipped=True, aligned=None, size_bytes=size)

    raise IngestError(
        f"Gzipped input {path} is neither BAM, FASTQ, nor FASTA (inner head: {inner_head[:8]!r})"
    )


def _bam_or_cram_has_alignment(path: Path) -> bool:
    """Return True if the first non-header record in the BAM/CRAM is mapped.

    Uses pysam if available; falls back to ``False`` on import failure so the
    function still works in environments that haven't installed htslib yet
    (e.g. early lint runs). Treating "unknown" as "unaligned" is conservative:
    the pipeline will then realign, which never gives a worse result.
    """
    try:
        import pysam
    except ImportError:
        return False
    try:
        with pysam.AlignmentFile(str(path), "rb", check_sq=False) as bam:
            for read in bam:
                return not read.is_unmapped
    except (ValueError, OSError):
        return False
    return False
