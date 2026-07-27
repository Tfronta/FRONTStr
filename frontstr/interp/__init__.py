"""Interpretation layer (Layer 3): forensic decisions.

Public surface:

- :func:`interpret_marker`, :func:`interpret_run` — orchestrators.
- :class:`Allele`, :class:`MarkerResult`, :class:`AlleleStatus`,
  :class:`CallRule`, :class:`TriType` — data model.
- :func:`build_expected_stutter`, :func:`classify_allele`,
  :func:`call_profile`, :func:`cross_check` — building blocks for tests
  and advanced callers.
- :func:`compress_isfg`, :func:`motif_repeat_summary`, :func:`ce_from_length`,
  :func:`ce_from_brackets` — nomenclature helpers.
"""

from frontstr.interp.classify import classify_allele
from frontstr.interp.concordance import cross_check
from frontstr.interp.isfg import (
    ce_from_brackets,
    ce_from_length,
    compress_isfg,
    motif_repeat_summary,
)
from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    MarkerResult,
    TriType,
)
from frontstr.interp.profile import (
    index_longtr_results,
    interpret_marker,
    interpret_run,
)
from frontstr.interp.stutter import (
    MotifRun,
    build_expected_stutter,
    find_motif_runs,
    top_two_runs,
    virtual_stutters,
)
from frontstr.interp.triallelic import call_profile

__all__ = [
    "Allele",
    "AlleleStatus",
    "CallRule",
    "MarkerResult",
    "MotifRun",
    "TriType",
    "build_expected_stutter",
    "call_profile",
    "ce_from_brackets",
    "ce_from_length",
    "classify_allele",
    "compress_isfg",
    "cross_check",
    "find_motif_runs",
    "index_longtr_results",
    "interpret_marker",
    "interpret_run",
    "motif_repeat_summary",
    "top_two_runs",
    "virtual_stutters",
]
