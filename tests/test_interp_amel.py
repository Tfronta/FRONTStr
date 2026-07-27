"""Tests for :mod:`frontstr.interp.amel` — AMEL sex typing."""

from __future__ import annotations

from pathlib import Path

import pysam

from frontstr.interp.amel import interpret_amel
from frontstr.interp.models import CallRule
from frontstr.panel.models import System

# Synthetic AMEL coordinates (0-based half-open in BAM; 1-based in System)
_CHROM_X = "chrX"
_CHROM_Y = "chrY"
_X_START = 1000  # ref_start (1-based) in System
_X_END = 1100  # ref_end (1-based) in System  → pileup window [999, 1100)
_Y_START = 2000
_Y_END = 2100
_CHR_LEN = 10_000
_FLANK = 50


def _amel_system(*, with_y: bool = True) -> System:
    s: dict[str, object] = {
        "name": "AMEL",
        "marker_type": "amel",
        "chromosome": _CHROM_X,
        "ref_start": _X_START,
        "ref_end": _X_END,
        "motif": "A",
        "period": 1,
        "category": "x_chromosomal",
    }
    if with_y:
        s["y_chromosome"] = _CHROM_Y
        s["y_ref_start"] = _Y_START
        s["y_ref_end"] = _Y_END
    return System(**s)


def _write_amel_bam(out_path: Path, x_reads: int, y_reads: int) -> Path:
    """Create a BAM with ``x_reads`` spanning [_X_START-1, _X_END) on chrX
    and ``y_reads`` spanning [_Y_START-1, _Y_END) on chrY."""
    header_dict = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": _CHROM_X, "LN": _CHR_LEN},
            {"SN": _CHROM_Y, "LN": _CHR_LEN},
        ],
        "RG": [{"ID": "rg1", "SM": "s", "LB": "l", "PL": "ONT"}],
    }
    raw = out_path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(raw), "wb", header=header_dict) as bam:
        for chrom, n_reads, ref_start_0 in [
            (_CHROM_X, x_reads, _X_START - 1),
            (_CHROM_Y, y_reads, _Y_START - 1),
        ]:
            ref_id = {_CHROM_X: 0, _CHROM_Y: 1}[chrom]
            tr_len = (_X_END - _X_START) if chrom == _CHROM_X else (_Y_END - _Y_START)
            seq_len = _FLANK + tr_len + _FLANK
            for i in range(n_reads):
                a = pysam.AlignedSegment(header=bam.header)
                a.query_name = f"{chrom[3]}{i}"
                a.query_sequence = "A" * seq_len
                a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
                a.reference_id = ref_id
                a.reference_start = ref_start_0 - _FLANK
                a.mapping_quality = 60
                a.flag = 0
                a.cigartuples = [(0, seq_len)]
                a.set_tag("RG", "rg1")
                bam.write(a)
    pysam.sort("-o", str(out_path), str(raw))
    pysam.index(str(out_path))
    raw.unlink(missing_ok=True)
    return out_path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_amel_female(tmp_path: Path) -> None:
    """Only X reads → female (X, HOMOZYGOUS)."""
    bam = _write_amel_bam(tmp_path / "female.bam", x_reads=20, y_reads=0)
    result = interpret_amel(_amel_system(), bam)
    assert result.call_rule == CallRule.HOMOZYGOUS
    assert len(result.alleles_called) == 1
    assert result.alleles_called[0].isfg == "X"
    assert result.alleles_called[0].ce is None


def test_amel_male(tmp_path: Path) -> None:
    """X and Y reads → male (X/Y, HETEROZYGOUS)."""
    bam = _write_amel_bam(tmp_path / "male.bam", x_reads=20, y_reads=15)
    result = interpret_amel(_amel_system(), bam)
    assert result.call_rule == CallRule.HETEROZYGOUS
    assert len(result.alleles_called) == 2
    isfg_calls = {a.isfg for a in result.alleles_called}
    assert isfg_calls == {"X", "Y"}
    assert all(a.ce is None for a in result.alleles_called)


def test_amel_no_reads(tmp_path: Path) -> None:
    """No reads → NO_DATA."""
    bam = _write_amel_bam(tmp_path / "empty.bam", x_reads=0, y_reads=0)
    result = interpret_amel(_amel_system(), bam)
    assert result.call_rule == CallRule.NO_DATA
    assert result.alleles_called == []
    assert result.total_reads == 0


def test_amel_below_min_reads_threshold(tmp_path: Path) -> None:
    """Y reads below min_reads (default 3) → treated as absent → female call."""
    bam = _write_amel_bam(tmp_path / "low_y.bam", x_reads=20, y_reads=2)
    result = interpret_amel(_amel_system(), bam)
    assert result.call_rule == CallRule.HOMOZYGOUS
    isfg_calls = {a.isfg for a in result.alleles_called}
    assert isfg_calls == {"X"}


def test_amel_without_y_region_config(tmp_path: Path) -> None:
    """System without Y coords → only chrX pileup → female regardless of BAM."""
    bam = _write_amel_bam(tmp_path / "no_y_cfg.bam", x_reads=20, y_reads=15)
    result = interpret_amel(_amel_system(with_y=False), bam)
    # Y never pileup'd → only X allele called
    assert result.call_rule == CallRule.HOMOZYGOUS
    assert result.alleles_called[0].isfg == "X"


def test_amel_total_reads_sums_both_chroms(tmp_path: Path) -> None:
    bam = _write_amel_bam(tmp_path / "total.bam", x_reads=20, y_reads=15)
    result = interpret_amel(_amel_system(), bam)
    assert result.total_reads == 35
