"""Export formats.

Currently implemented (Phase 1.7):

- :func:`write_run_json` — full machine-readable payload as JSON.
- :func:`write_profile_csv` — wide-format forensic profile (1 row x marker).
- :func:`write_evidence_csv` — long-format cluster evidence (1 row x cluster).
- :func:`write_seqs_csv` — per-allele ISFG + consensus trail.
- :func:`write_run_vcf` — native, sequence-resolved VCF (ALT is the allele
  sequence, so iso-alleles survive the round trip).
- :func:`write_run_xlsx` — multi-sheet review workbook.
- :func:`write_tidy` — cohort-scale long table (one row per sample x marker x
  allele) as CSV and Parquet, built from run JSONs.

Pending: CODIS CMF, NIST MIDST, PDF, ZIP bundle.
"""

from frontstr.exports.csv import (
    EVIDENCE_HEADERS,
    PROFILE_HEADERS,
    SEQS_HEADERS,
    write_evidence_csv,
    write_profile_csv,
    write_seqs_csv,
)
from frontstr.exports.json import JsonMode, write_run_json
from frontstr.exports.tidy import build_tidy_rows, load_payloads, write_tidy
from frontstr.exports.vcf import write_run_vcf
from frontstr.exports.xlsx import write_run_xlsx

__all__ = [
    "EVIDENCE_HEADERS",
    "PROFILE_HEADERS",
    "SEQS_HEADERS",
    "JsonMode",
    "build_tidy_rows",
    "load_payloads",
    "write_evidence_csv",
    "write_profile_csv",
    "write_run_json",
    "write_run_vcf",
    "write_run_xlsx",
    "write_seqs_csv",
    "write_tidy",
]
