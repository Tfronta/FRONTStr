"""Tests for :func:`frontstr.evidence.pileup.pileup_locus`."""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from frontstr.errors import EvidenceError
from frontstr.evidence.pileup import pileup_locus
from tests.conftest import SYNTH_CHR_LEN, SYNTH_CHROM, SYNTH_FLANK_LEN, SYNTH_TR_END, SYNTH_TR_START

# ---------------------------------------------------------------------------
# Helpers for CRAM tests
# ---------------------------------------------------------------------------


def _write_synth_fasta(path: Path) -> Path:
    """Write a minimal FASTA for SYNTH_CHROM and index it with pysam.faidx."""
    seq = "A" * SYNTH_CHR_LEN
    path.write_text(f">{SYNTH_CHROM}\n{seq}\n")
    pysam.faidx(str(path))
    return path


def _bam_to_cram(bam_path: Path, cram_path: Path, fasta_path: Path) -> Path:
    """Convert an existing sorted+indexed BAM to CRAM using the given reference."""
    pysam.view(
        "-C",
        "-T",
        str(fasta_path),
        "-o",
        str(cram_path),
        str(bam_path),
        catch_stdout=False,
    )
    pysam.index(str(cram_path))
    return cram_path


def _make_boundary_bam(out_path: Path, cigartuples: list[tuple[int, int]], seq: str) -> Path:
    """Write a single-read BAM with the given CIGAR, starting at ref pos 0."""
    header_dict = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": SYNTH_CHROM, "LN": SYNTH_CHR_LEN}],
        "RG": [{"ID": "rg1", "SM": "s", "LB": "l", "PL": "ONT"}],
    }
    raw = out_path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(raw), "wb", header=header_dict) as bam:
        a = pysam.AlignedSegment(header=bam.header)
        a.query_name = "del_read"
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.flag = 0
        a.cigartuples = cigartuples
        a.set_tag("RG", "rg1")
        bam.write(a)
    pysam.sort("-o", str(out_path), str(raw))
    pysam.index(str(out_path))
    raw.unlink(missing_ok=True)
    return out_path


def test_pileup_heterozygous_counts(synth_bam_heterozygous: Path) -> None:
    obs = pileup_locus(synth_bam_heterozygous, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    assert len(obs) == 11

    lengths = sorted(len(o.sequence) for o in obs)
    assert lengths.count(48) == 5  # CE 12
    assert lengths.count(44) == 4  # CE 11
    assert lengths.count(40) == 2  # CE 10


def test_pileup_haplotype_tags(synth_bam_heterozygous: Path) -> None:
    obs = pileup_locus(synth_bam_heterozygous, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    hp1 = [o for o in obs if o.hp == 1]
    hp2 = [o for o in obs if o.hp == 2]
    hp_none = [o for o in obs if o.hp is None]
    assert len(hp1) == 5
    assert len(hp2) == 4
    assert len(hp_none) == 2

    assert all(len(o.sequence) == 48 for o in hp1)
    assert all(len(o.sequence) == 44 for o in hp2)
    assert all(len(o.sequence) == 40 for o in hp_none)


def test_pileup_sequences_are_motif_repeats(synth_bam_heterozygous: Path) -> None:
    obs = pileup_locus(synth_bam_heterozygous, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    for o in obs:
        n = len(o.sequence) // 4
        assert o.sequence == "AGAT" * n, f"unexpected sequence for {o.read_id}"


def test_pileup_mean_quality_set(synth_bam_heterozygous: Path) -> None:
    obs = pileup_locus(synth_bam_heterozygous, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    for o in obs:
        # 'I' ASCII = 73; Phred = 73 - 33 = 40
        assert 39.0 <= o.mean_qual <= 40.0


def test_pileup_drops_low_mapq(synth_bam_low_mapq: Path) -> None:
    obs = pileup_locus(synth_bam_low_mapq, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END, min_mapq=20)
    assert len(obs) == 3
    assert {o.read_id for o in obs} == {"hi0", "hi1", "hi2"}


def test_pileup_keeps_low_mapq_when_threshold_lowered(synth_bam_low_mapq: Path) -> None:
    obs = pileup_locus(synth_bam_low_mapq, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END, min_mapq=0)
    assert len(obs) == 6


def test_pileup_strand_tracking(synth_bam_homozygous: Path) -> None:
    obs = pileup_locus(synth_bam_homozygous, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    forward = [o for o in obs if o.strand == "+"]
    reverse = [o for o in obs if o.strand == "-"]
    assert len(forward) == 4
    assert len(reverse) == 4
    for o in obs:
        assert o.sequence == "AGAT" * 12


def test_pileup_rejects_invalid_window(tmp_path: Path) -> None:
    bam = tmp_path / "x.bam"
    bam.write_bytes(b"")
    with pytest.raises(EvidenceError, match="not found"):
        pileup_locus(tmp_path / "missing.bam", "chr1", 100, 200)
    with pytest.raises(EvidenceError, match="Invalid locus window"):
        pileup_locus(tmp_path / "missing.bam", "chr1", 100, 100)
    with pytest.raises(EvidenceError, match="Invalid locus window"):
        pileup_locus(tmp_path / "missing.bam", "chr1", -1, 100)


def test_pileup_missing_bam(tmp_path: Path) -> None:
    with pytest.raises(EvidenceError, match="not found"):
        pileup_locus(tmp_path / "nope.bam", "chr1", 100, 200)


# ---------------------------------------------------------------------------
# Bug #3 regression tests: deletion exactly at the TR boundary
# ---------------------------------------------------------------------------


def test_deletion_at_tr_start_not_dropped(tmp_path: Path) -> None:
    """A read with a deletion at rpos==start must not be silently dropped.

    CIGAR: 100M 4D 44M 100M
    The first 4 bp of the TR (ref 100-103) are deleted in this read.
    The read should be returned with a 44-bp sequence (11 AGAT repeats).
    """
    del_bp = 4
    kept_tr = SYNTH_TR_END - SYNTH_TR_START - del_bp  # 44
    seq = "A" * SYNTH_FLANK_LEN + "AGAT" * (kept_tr // 4) + "T" * SYNTH_FLANK_LEN
    cigar = [
        (0, SYNTH_FLANK_LEN),  # 100M left flank
        (2, del_bp),  # 4D — deletion at rpos == SYNTH_TR_START
        (0, kept_tr),  # 44M — rest of TR
        (0, SYNTH_FLANK_LEN),  # 100M right flank
    ]
    bam = _make_boundary_bam(tmp_path / "del_start.bam", cigar, seq)
    obs = pileup_locus(bam, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    assert len(obs) == 1, "read with deletion at TR start must not be dropped"
    assert len(obs[0].sequence) == kept_tr


def test_deletion_at_tr_end_not_dropped(tmp_path: Path) -> None:
    """A read with a deletion at rpos==end must not be silently dropped.

    CIGAR: 100M 48M 4D 96M
    Ref positions 148-151 are deleted (148 == SYNTH_TR_END, the exclusive end).
    The full TR (48 bp) is present in the read; the deletion is in the right flank.
    """
    del_bp = 4
    full_tr = SYNTH_TR_END - SYNTH_TR_START  # 48
    right_flank = SYNTH_FLANK_LEN - del_bp  # 96
    seq = "A" * SYNTH_FLANK_LEN + "AGAT" * (full_tr // 4) + "T" * right_flank
    cigar = [
        (0, SYNTH_FLANK_LEN),  # 100M left flank
        (0, full_tr),  # 48M full TR
        (2, del_bp),  # 4D — deletion at rpos == SYNTH_TR_END
        (0, right_flank),  # 96M right flank (shortened by deletion)
    ]
    bam = _make_boundary_bam(tmp_path / "del_end.bam", cigar, seq)
    obs = pileup_locus(bam, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    assert len(obs) == 1, "read with deletion at TR end must not be dropped"
    assert len(obs[0].sequence) == full_tr


# ---------------------------------------------------------------------------
# CRAM support tests
# ---------------------------------------------------------------------------


def test_cram_requires_reference_fasta(tmp_path: Path) -> None:
    """pileup_locus on a .cram file without reference_fasta raises EvidenceError."""
    cram = tmp_path / "dummy.cram"
    cram.write_bytes(b"")
    with pytest.raises(EvidenceError, match="reference FASTA"):
        pileup_locus(cram, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)


def test_cram_pileup_returns_same_observations(
    synth_bam_heterozygous: Path, tmp_path: Path
) -> None:
    """A CRAM converted from the synthetic BAM must yield identical observations."""
    fasta = _write_synth_fasta(tmp_path / "ref.fa")
    cram = _bam_to_cram(synth_bam_heterozygous, tmp_path / "het.cram", fasta)

    obs_bam = pileup_locus(synth_bam_heterozygous, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END)
    obs_cram = pileup_locus(cram, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END, reference_fasta=fasta)

    assert len(obs_cram) == len(obs_bam)
    bam_lengths = sorted(len(o.sequence) for o in obs_bam)
    cram_lengths = sorted(len(o.sequence) for o in obs_cram)
    assert cram_lengths == bam_lengths
