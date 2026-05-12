"""BAM/CRAM header and content validation prior to handing off to LongTR.

LongTR is strict about read groups (it uses ``@RG SM``/``LB``) and is happy to
silently produce nonsense if the reference build mismatches the alignment.
This module enforces a small set of checks and returns a list of human-readable
warnings; fatal issues raise :class:`IngestError`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from frontstr.errors import IngestError


@dataclass(slots=True)
class ValidationReport:
    """Outcome of :func:`validate_bam`.

    ``warnings`` is a non-fatal list (the pipeline can still proceed but
    we surface these in the run's ``ingest_warnings`` audit field).
    """

    warnings: list[str] = field(default_factory=list)
    n_records_inspected: int = 0
    median_mapq: int = 0
    sample_names: list[str] = field(default_factory=list)
    reference_names: list[str] = field(default_factory=list)


_MAPQ_SAMPLE_SIZE = 10_000
_LOW_MAPQ_THRESHOLD = 10


def validate_bam(path: Path, expected_build: str | None = None) -> ValidationReport:
    """Validate a BAM/CRAM that is being passed through to LongTR.

    Args:
        path: Path to an indexed (or indexable) BAM/CRAM.
        expected_build: e.g. ``"GRCh38"``; used to compare ``@SQ`` M5 hashes
            when a reference MD5 table is available. ``None`` skips that check.

    Returns:
        :class:`ValidationReport` with collected warnings and stats.

    Raises:
        IngestError: For fatal issues (missing ``@RG``, missing ``SM``).
    """
    try:
        import pysam
    except ImportError as exc:
        raise IngestError(
            "pysam is required for BAM validation; install with `pip install pysam`"
        ) from exc

    report = ValidationReport()

    try:
        bam = pysam.AlignmentFile(str(path), "rb")
    except (ValueError, OSError) as exc:
        raise IngestError(f"Cannot open BAM/CRAM {path}: {exc}") from exc

    try:
        header = bam.header.to_dict()
    except Exception as exc:
        raise IngestError(f"Cannot parse header of {path}: {exc}") from exc

    rgs = header.get("RG", [])
    if not rgs:
        raise IngestError(
            f"BAM {path} is missing @RG header. "
            "Use `samtools addreplacerg` to add a read group with SM/LB."
        )
    missing_sm = [rg for rg in rgs if "SM" not in rg]
    if missing_sm:
        raise IngestError(
            f"BAM {path} has @RG entries without SM tag: {missing_sm!r}. "
            "LongTR requires SM to assign reads to samples."
        )
    report.sample_names = sorted({rg["SM"] for rg in rgs})

    sqs = header.get("SQ", [])
    report.reference_names = [sq["SN"] for sq in sqs]
    if not sqs:
        report.warnings.append("@SQ section is empty; reference is unknown")

    if expected_build:
        report.warnings.extend(_validate_reference_md5(sqs, expected_build))

    bai = path.with_suffix(path.suffix + ".bai")
    csi = path.with_suffix(path.suffix + ".csi")
    if not bai.exists() and not csi.exists():
        report.warnings.append("BAM index missing; will attempt to create one")
        try:
            pysam.index(str(path))
        except pysam.SamtoolsError as exc:
            raise IngestError(f"Cannot index BAM {path}: {exc}") from exc

    mapqs: list[int] = []
    for i, read in enumerate(bam.fetch(until_eof=True)):
        if i >= _MAPQ_SAMPLE_SIZE:
            break
        if not read.is_unmapped:
            mapqs.append(read.mapping_quality)
    report.n_records_inspected = len(mapqs)
    report.median_mapq = _median(mapqs) if mapqs else 0
    if mapqs and report.median_mapq < _LOW_MAPQ_THRESHOLD:
        report.warnings.append(
            f"Median MAPQ={report.median_mapq} (< {_LOW_MAPQ_THRESHOLD}); "
            "LongTR may drop many reads"
        )

    bam.close()
    return report


def _validate_reference_md5(sqs: list[dict], expected_build: str) -> list[str]:
    """Compare @SQ M5 fields with a known reference MD5 table.

    The actual MD5 table is loaded lazily from :mod:`frontstr.panel.reference`
    once implemented; until then we return an empty list (no false positives).
    """
    # TODO(phase-1.2): wire in frontstr.panel.reference.MD5_TABLES[build].
    _ = sqs, expected_build
    return []


def _median(values: list[int]) -> int:
    """Median of a non-empty list of ints, biased to the lower middle."""
    if not values:
        return 0
    s = sorted(values)
    return s[len(s) // 2]
