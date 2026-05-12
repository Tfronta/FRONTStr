"""LongTR caller — argv builder, subprocess wrapper, high-level orchestration.

This is Layer 1 of FRONTStr. It exists to (a) call LongTR with the right
arguments for the platform and (b) hand its VCF off to the parser. It is
deliberately thin: forensic decisions never live here.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from frontstr.caller.bed import write_panel_bed
from frontstr.caller.vcf import LongTRResult, parse_longtr_vcf
from frontstr.errors import CallerError
from frontstr.panel.models import Panel

_DEFAULT_BINARY = "LongTR"
_BINARY_ENV = "FRONTSTR_LONGTR_BIN"


# ---------------------------------------------------------------------------
# Argv builder
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class LongTRInvocation:
    """Resolved arguments and paths for a LongTR call."""

    argv: list[str]
    bam: Path
    bed: Path
    vcf: Path


def build_longtr_argv(
    *,
    bam: Path,
    fasta: Path,
    bed: Path,
    chrom: str | None,
    vcf_out: Path,
    platform: str,
    min_reads: int = 10,
    max_tr_len: int = 1000,
    min_mapq: int = 20,
    indel_flank_len: int = 5,
    phased: bool = False,
    skip_assembly: bool = False,
    binary: str | None = None,
    extra: list[str] | None = None,
) -> LongTRInvocation:
    """Construct a deterministic LongTR argv for the given inputs.

    See plan-longtr-improved.md §4 and §18 for parameter rationale.

    Args:
        bam: Indexed BAM (single-sample for forensic use).
        fasta: Indexed reference FASTA.
        bed: BED file of TR regions (see :func:`frontstr.caller.bed.write_panel_bed`).
        chrom: Restrict to one chromosome via ``--chrom``; ``None`` to process all.
        vcf_out: Output ``.vcf.gz`` path.
        platform: ``"ont"`` (R10 simplex defaults) or ``"hifi"``.
        binary: Override the ``LongTR`` binary location.
        extra: Extra args appended verbatim (advanced use).
    """
    if platform not in {"ont", "hifi"}:
        raise CallerError(f"Unsupported platform: {platform!r}")

    argv: list[str] = [
        binary or _resolve_binary(),
        "--bams", str(bam),
        "--fasta", str(fasta),
        "--regions", str(bed),
        "--tr-vcf", str(vcf_out),
        "--min-reads", str(min_reads),
        "--max-tr-len", str(max_tr_len),
        "--min-mapq", str(min_mapq),
        "--indel-flank-len", str(indel_flank_len),
        "--lib-from-samp",
    ]
    if chrom:
        argv += ["--chrom", chrom]
    if platform == "ont":
        argv += [
            "--min-mean-qual", "10",
            "--alignment-params",
            "-1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1.0,-1.0",
        ]
    else:
        argv += ["--min-mean-qual", "30"]
    if phased:
        argv += ["--phased-bam"]
    if skip_assembly:
        argv += ["--skip-assembly"]
    if extra:
        argv += list(extra)
    return LongTRInvocation(argv=argv, bam=bam, bed=bed, vcf=vcf_out)


def _resolve_binary() -> str:
    """Find the LongTR binary or raise."""
    binary = os.environ.get(_BINARY_ENV) or _DEFAULT_BINARY
    if shutil.which(binary) is None:
        raise CallerError(
            f"LongTR binary {binary!r} not found on PATH. "
            f"Set {_BINARY_ENV} or install LongTR (https://github.com/gymrek-lab/LongTR)."
        )
    return binary


# ---------------------------------------------------------------------------
# Subprocess runner
# ---------------------------------------------------------------------------


def run_longtr(
    inv: LongTRInvocation,
    *,
    log_path: Path | None = None,
    timeout_s: float | None = None,
) -> Path:
    """Execute LongTR. Returns the VCF path on success.

    Args:
        inv: Output of :func:`build_longtr_argv`.
        log_path: Optional file to capture combined stdout/stderr.
        timeout_s: Wall-clock timeout in seconds.

    Raises:
        CallerError: On non-zero exit or timeout. The exception message
            includes the tail of the captured log when available.
    """
    inv.vcf.parent.mkdir(parents=True, exist_ok=True)
    log_fh: Any
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_fh = log_path.open("w", encoding="utf-8")
    else:
        log_fh = subprocess.PIPE

    try:
        proc = subprocess.run(
            inv.argv,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=timeout_s,
            text=True,
        )
    except subprocess.TimeoutExpired as exc:
        raise CallerError(
            f"LongTR timed out after {timeout_s}s: {' '.join(inv.argv[:3])} ..."
        ) from exc
    except FileNotFoundError as exc:
        raise CallerError(f"LongTR binary not found: {inv.argv[0]}") from exc
    finally:
        if log_path is not None:
            log_fh.close()

    if proc.returncode != 0:
        tail = _last_lines(log_path) if log_path else (proc.stdout or "")[-2000:]
        raise CallerError(
            f"LongTR exited with code {proc.returncode}. Last output:\n{tail.strip()}"
        )
    if not inv.vcf.exists():
        raise CallerError(f"LongTR finished but did not write {inv.vcf}")
    return inv.vcf


def _last_lines(path: Path | None, n: int = 40) -> str:
    if path is None or not path.exists():
        return ""
    lines = path.read_text(errors="replace").splitlines()
    return "\n".join(lines[-n:])


# ---------------------------------------------------------------------------
# High-level orchestrator
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class LongTRRun:
    """One LongTR execution: paths produced + parsed results."""

    bed: Path
    vcf: Path
    log: Path | None
    results: list[LongTRResult]


class LongTRRunner:
    """Run LongTR end-to-end for a panel and parse its VCF.

    Typical use::

        runner = LongTRRunner(panel=panel, reference=ref, platform="ont")
        run = runner.run(bam=bam_path, out_dir=outdir)
        for r in run.results:
            print(r.marker_name, r.alleles)
    """

    def __init__(
        self,
        *,
        panel: Panel,
        reference: Path,
        platform: str,
        min_reads: int = 10,
        max_tr_len: int = 1000,
        min_mapq: int = 20,
        indel_flank_len: int = 5,
        phased: bool = False,
        skip_assembly: bool = False,
        binary: str | None = None,
    ) -> None:
        self.panel = panel
        self.reference = reference
        self.platform = platform
        self.min_reads = min_reads
        self.max_tr_len = max_tr_len
        self.min_mapq = min_mapq
        self.indel_flank_len = indel_flank_len
        self.phased = phased
        self.skip_assembly = skip_assembly
        self.binary = binary

    def run(
        self,
        *,
        bam: Path,
        out_dir: Path,
        chrom: str | None = None,
        timeout_s: float | None = None,
    ) -> LongTRRun:
        """Write BED, invoke LongTR, parse the VCF.

        Args:
            bam: Indexed sample BAM.
            out_dir: Where to write the BED, VCF, and log.
            chrom: Restrict to a single chromosome (default: all).
            timeout_s: LongTR wall-clock timeout.
        """
        out_dir.mkdir(parents=True, exist_ok=True)
        bed = write_panel_bed(self.panel, out_dir / "panel.bed")
        vcf = out_dir / "tr_calls.vcf.gz"
        log = out_dir / "longtr.log"

        inv = build_longtr_argv(
            bam=bam,
            fasta=self.reference,
            bed=bed,
            chrom=chrom,
            vcf_out=vcf,
            platform=self.platform,
            min_reads=self.min_reads,
            max_tr_len=self.max_tr_len,
            min_mapq=self.min_mapq,
            indel_flank_len=self.indel_flank_len,
            phased=self.phased,
            skip_assembly=self.skip_assembly,
            binary=self.binary,
        )
        run_longtr(inv, log_path=log, timeout_s=timeout_s)
        results = parse_longtr_vcf(vcf)
        return LongTRRun(bed=bed, vcf=vcf, log=log, results=results)
