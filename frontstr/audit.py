"""Audit trail: what was run, with what, and what it flagged.

A forensic result is only as good as the account of how it was produced. Four
layers, in increasing durability:

0. **Locus trace** — :mod:`frontstr.trace`, written per sample as
   ``<sample>.trace.txt``. Answers the question the other three cannot:
   *where* did a call go wrong. Every read rejected and why, the length bins,
   the clusters and their consensus, the aligned sequences, the haplotype
   counts per candidate, how each allele number was derived, and why each
   candidate was called or discarded. It is the difference between a caller
   that can be audited and a black box, so ``frontstr batch`` writes it by
   default rather than on request.
1. **Process log** — a JSONL record of what the pipeline did, one event per
   line, written next to the outputs. Answers "what happened, in what order,
   and how long did it take". Ephemeral in the sense that it is not part of the
   result, but it is what you read when a run behaves oddly.
2. **Audit record** — :class:`AuditRecord`, embedded in the canonical JSON.
   Answers "under exactly what configuration was this result produced": tool
   version, POA backend, stutter model, QC thresholds, input hashes, and the
   flag census. This travels *with* the result and is what a reviewer or an
   opposing expert reads.
3. **Integrity hash** — a SHA-256 over the audit record's own canonical form,
   so a record cannot be edited without detection.

Why the flag census lives here
------------------------------

A reviewer needs to know not just that a marker was flagged, but that the run
as a whole raised (say) four warnings and no errors — and, just as importantly,
*which conditions were checked at all*. :attr:`AuditRecord.flags_checked` lists
every code the pipeline is capable of raising, so a code absent from
:attr:`AuditRecord.flag_counts` provably means "checked and not found" rather
than "never looked at". That distinction is the difference between a clean
report and an incomplete one.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

from pydantic import BaseModel, Field

from frontstr.interp.models import FlagCode, MarkerResult
from frontstr.interp.qc import QcThresholds
from frontstr.log import PROCESS_LOG_NAME, configure_logging, get_logger
from frontstr.version import __version__

#: Bump when the audit record's shape changes incompatibly.
AUDIT_SCHEMA_VERSION = "1.0"

__all__ = [
    "AUDIT_SCHEMA_VERSION",
    "PROCESS_LOG_NAME",
    "AuditRecord",
    "InputFile",
    "build_audit_record",
    "configure_logging",
    "file_sha256",
    "get_logger",
]


# ---------------------------------------------------------------------------
# Layer 2 — audit record
# ---------------------------------------------------------------------------


class InputFile(BaseModel):
    """One input, with the hash it had when it was read."""

    role: str
    path: str
    sha256: str | None = None


class AuditRecord(BaseModel):
    """Everything needed to reproduce and to challenge a result."""

    schema_version: str = AUDIT_SCHEMA_VERSION
    tool_name: str = "frontstr"
    tool_version: str = __version__

    #: Resolved software the numbers actually depend on — a result produced
    #: without a POA backend is not the same result.
    poa_backend: str = ""
    stutter_model_version: str = ""
    stutter_model_protocol: str = ""

    #: Every threshold that moved a call, in one place.
    analytical_thresh: float | None = None
    calling_thresh: float | None = None
    qc_thresholds: QcThresholds = Field(default_factory=QcThresholds)

    inputs: list[InputFile] = Field(default_factory=list)

    #: Flag codes raised in this run, with counts.
    flag_counts: dict[str, int] = Field(default_factory=dict)
    #: Counts by severity, for the one-line verdict.
    severity_counts: dict[str, int] = Field(default_factory=dict)
    #: Every code the pipeline can raise. A code listed here but absent from
    #: ``flag_counts`` was checked and not found — the guarantee that makes a
    #: clean report meaningful.
    flags_checked: list[str] = Field(default_factory=list)
    #: Markers carrying at least one WARN or ERROR, for triage.
    markers_needing_review: list[str] = Field(default_factory=list)

    #: SHA-256 over this record's canonical JSON, excluding this field.
    integrity_sha256: str | None = None

    def sealed(self) -> AuditRecord:
        """Return a copy carrying its own integrity hash.

        The hash covers the canonical JSON of every other field, so editing any
        of them after the fact no longer matches. It proves internal
        consistency, not authenticity — it is a tamper-evidence measure, not a
        signature, and nothing here pretends otherwise.
        """
        body = self.model_dump(mode="json", exclude={"integrity_sha256"})
        canonical = json.dumps(body, sort_keys=True, separators=(",", ":"))
        digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
        return self.model_copy(update={"integrity_sha256": digest})

    def verify(self) -> bool:
        """True when :attr:`integrity_sha256` still matches the record body."""
        if self.integrity_sha256 is None:
            return False
        return self.sealed().integrity_sha256 == self.integrity_sha256


def build_audit_record(
    results: list[MarkerResult],
    *,
    inputs: list[InputFile] | None = None,
    qc_thresholds: QcThresholds | None = None,
    analytical_thresh: float | None = None,
    calling_thresh: float | None = None,
) -> AuditRecord:
    """Assemble the audit record for a finished run and seal it.

    Resolves the POA backend and stutter model at call time rather than taking
    them as arguments: the point is to record what the run *actually used*, and
    a caller passing those in could pass something else.
    """
    from frontstr.evidence.consensus import poa_backend_name
    from frontstr.panel.stutter_calib import DEFAULT_STUTTER_MODEL

    flag_counts: Counter[str] = Counter()
    severity_counts: Counter[str] = Counter()
    needs_review: list[str] = []
    for r in results:
        flagged = False
        for f in r.flags:
            flag_counts[f.code.value] += 1
            severity_counts[f.severity.value] += 1
            if f.severity.value in ("warn", "error"):
                flagged = True
        for a in r.alleles:
            for f in a.flags:
                flag_counts[f.code.value] += 1
                severity_counts[f.severity.value] += 1
                if f.severity.value in ("warn", "error"):
                    flagged = True
        if flagged:
            needs_review.append(r.marker_name)

    record = AuditRecord(
        poa_backend=poa_backend_name() or "none",
        stutter_model_version=DEFAULT_STUTTER_MODEL.version,
        stutter_model_protocol=DEFAULT_STUTTER_MODEL.protocol,
        analytical_thresh=analytical_thresh,
        calling_thresh=calling_thresh,
        qc_thresholds=qc_thresholds or QcThresholds(),
        inputs=inputs or [],
        flag_counts=dict(sorted(flag_counts.items())),
        severity_counts=dict(sorted(severity_counts.items())),
        flags_checked=sorted(c.value for c in FlagCode),
        markers_needing_review=needs_review,
    )
    return record.sealed()


def file_sha256(path: Path, chunk: int = 1 << 20) -> str:
    """Stream a file through SHA-256 (1 MiB chunks)."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for buf in iter(lambda: fh.read(chunk), b""):
            h.update(buf)
    return h.hexdigest()
