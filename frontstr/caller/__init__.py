"""LongTR caller integration (Layer 1).

Public surface:

- :class:`LongTRRunner` — full orchestration (BED → invoke → parse).
- :func:`build_longtr_argv` / :func:`run_longtr` — low-level building blocks.
- :func:`parse_longtr_vcf` — standalone parser for an existing VCF.
- :class:`LongTRResult`, :class:`LongTRAlleleSpec`, :class:`LongTRSampleCall`.
- :func:`write_panel_bed`, :func:`split_panel_by_chromosome`.
"""

from frontstr.caller.longtr import (
    LongTRInvocation,
    LongTRRun,
    LongTRRunner,
    build_longtr_argv,
    run_longtr,
)
from frontstr.caller.vcf import (
    LongTRAlleleSpec,
    LongTRResult,
    LongTRSampleCall,
    parse_longtr_vcf,
)
from frontstr.panel.bed import split_panel_by_chromosome, write_panel_bed

__all__ = [
    "LongTRAlleleSpec",
    "LongTRInvocation",
    "LongTRResult",
    "LongTRRun",
    "LongTRRunner",
    "LongTRSampleCall",
    "build_longtr_argv",
    "parse_longtr_vcf",
    "run_longtr",
    "split_panel_by_chromosome",
    "write_panel_bed",
]
