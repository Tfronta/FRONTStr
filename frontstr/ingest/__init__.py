"""Input ingestion: format detection, BAM validation, alignment.

Public surface:

- :func:`detect_input` — magic-byte sniff (no extension dependency).
- :func:`validate_bam` — header / index / MAPQ sanity for pre-aligned BAM.
- :func:`align_to_reference` — minimap2 + samtools sort+index for FASTQ inputs.
"""

from frontstr.ingest.detect import InputKind, detect_input
from frontstr.ingest.validate import validate_bam

__all__ = ["InputKind", "detect_input", "validate_bam"]
