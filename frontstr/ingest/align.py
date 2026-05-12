"""Alignment with minimap2 + samtools sort/index.

Implemented in Phase 1.1 of ROADMAP.md.
"""

from __future__ import annotations

from pathlib import Path

from frontstr.errors import IngestError


def align_to_reference(
    input_path: Path,
    reference: Path,
    out_bam: Path,
    *,
    platform: str = "ont",
    threads: int = 4,
    sample: str = "S",
    library: str | None = None,
) -> Path:
    """Align FASTQ/uBAM input to ``reference`` and produce a sorted, indexed BAM.

    Args:
        input_path: FASTQ(.gz) or uBAM.
        reference: Indexed FASTA reference.
        out_bam: Destination BAM path.
        platform: ``"ont"`` (map-ont) or ``"hifi"`` (map-hifi).
        threads: Threads for minimap2 and samtools sort.
        sample: ``SM`` value injected into ``@RG``.
        library: ``LB`` value; defaults to ``sample`` if not provided.

    Returns:
        Path to the produced BAM (also writes ``.bai`` index).
    """
    _ = input_path, reference, out_bam, platform, threads, sample, library
    raise IngestError("align_to_reference is not yet implemented (Phase 1.1)")
