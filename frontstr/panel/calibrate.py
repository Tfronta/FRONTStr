"""Panel calibration: compute corr_value for each marker from a reference FASTA.

For period > 0 markers (ce_from_length):
    corr_value = window_bp - (longest_ref_repeat_run_bp)
    Analogue of Método B: flanking bp = window minus the repeat array.

For period = -1 markers (ce_from_brackets):
    corr_value = motif units in the window that fall OUTSIDE the core repeat block.
    The core block is the maximal group of motif occurrences where consecutive
    tokens are separated by <= GAP_THRESHOLD non-motif characters.  Legitimate
    intra-allele spacers (ta, tca, tccata, tcg, tccag …) are all <= 6 chars,
    so GAP_THRESHOLD=12 cleanly separates core from flanking noise.

    Calibration is done against the GRCh38 reference sequence (not sample reads)
    so the result is independent of any particular sample.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from frontstr.interp.isfg import ce_from_brackets, compress_isfg
from frontstr.panel.models import Panel, System

_GAP_THRESHOLD = 12  # chars of non-motif text; above this = flank boundary
_MOTIF_TOKEN_RE = re.compile(r"(\[[A-Z]+\]\d+|[A-Z]{2,})")


@dataclass
class CalibrationResult:
    marker_name: str
    period: int
    window_bp: int
    corr_value: int
    # period=-1 detail
    ref_isfg: str = ""
    raw_ce: float | None = None
    core_units: int | None = None
    flanking_units: int | None = None
    note: str = ""


def _rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _fetch_ref(fasta_path: Path, chrom: str, start0: int, end0: int) -> str:
    """Return uppercase reference sequence for [start0, end0) (0-based half-open)."""
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required for calibration") from exc
    with pysam.FastaFile(str(fasta_path)) as fa:
        return fa.fetch(chrom, start0, end0).upper()


def _motif_tokens(isfg: str) -> list[tuple[int, int, int]]:
    """Return (start, end, unit_count) for every motif token in isfg."""
    out = []
    for m in _MOTIF_TOKEN_RE.finditer(isfg):
        s = m.group()
        count_match = re.search(r"\d+", s) if s.startswith("[") else None
        n = int(count_match.group()) if count_match else 1
        out.append((m.start(), m.end(), n))
    return out


def _core_and_flanking(isfg: str, gap_threshold: int = _GAP_THRESHOLD) -> tuple[int, int]:
    """Return (core_units, flanking_units) for a period=-1 ISFG string.

    Groups motif tokens separated by <= gap_threshold chars into a single run.
    The run with the most units is the core repeat block; the remaining units
    are flanking noise caused by chance motif matches in the window flanks.
    """
    tokens = _motif_tokens(isfg)
    if not tokens:
        return 0, 0

    groups: list[list[tuple[int, int, int]]] = []
    cur: list[tuple[int, int, int]] = [tokens[0]]
    for tok in tokens[1:]:
        gap = tok[0] - cur[-1][1]
        if gap <= gap_threshold:
            cur.append(tok)
        else:
            groups.append(cur)
            cur = [tok]
    groups.append(cur)

    def _units(g: list[tuple[int, int, int]]) -> int:
        return sum(t[2] for t in g)

    core_units = max(_units(g) for g in groups)
    total_units = sum(_units(g) for g in groups)
    return core_units, total_units - core_units


def calibrate_marker_neg1(system: System, fasta_path: Path) -> CalibrationResult:
    """Calibrate corr_value for a period=-1 marker against the reference FASTA.

    corr_value = number of motif units in the extraction window that lie
    outside the main repeat block, as determined from the GRCh38 reference
    sequence.  These flanking occurrences inflate ce_from_brackets and must
    be subtracted.
    """
    if system.period != -1:
        raise ValueError(f"{system.name}: expected period=-1, got {system.period}")

    start0 = system.ref_start - 1
    end0 = system.ref_end
    window_bp = end0 - start0

    ref_seq = _fetch_ref(fasta_path, system.chromosome, start0, end0)
    ref_isfg = compress_isfg(ref_seq, motif=system.motif, strand=system.strand)
    raw_ce = ce_from_brackets(ref_isfg)
    core_units, flanking_units = _core_and_flanking(ref_isfg)

    return CalibrationResult(
        marker_name=system.name,
        period=-1,
        window_bp=window_bp,
        corr_value=flanking_units,
        ref_isfg=ref_isfg,
        raw_ce=raw_ce,
        core_units=core_units,
        flanking_units=flanking_units,
    )


def calibrate_panel_neg1(
    panel: Panel,
    fasta_path: Path,
) -> list[CalibrationResult]:
    """Calibrate all period=-1 STR markers in *panel* against the reference FASTA."""
    return [
        calibrate_marker_neg1(s, fasta_path)
        for s in panel.systems
        if s.period == -1 and s.marker_type == "str"
    ]
