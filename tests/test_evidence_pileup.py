"""Tests for :func:`frontstr.evidence.pileup.pileup_locus`."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.errors import EvidenceError
from frontstr.evidence.pileup import pileup_locus
from tests.conftest import SYNTH_CHROM, SYNTH_TR_END, SYNTH_TR_START


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
    obs = pileup_locus(
        synth_bam_low_mapq, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END, min_mapq=20
    )
    assert len(obs) == 3
    assert {o.read_id for o in obs} == {"hi0", "hi1", "hi2"}


def test_pileup_keeps_low_mapq_when_threshold_lowered(synth_bam_low_mapq: Path) -> None:
    obs = pileup_locus(
        synth_bam_low_mapq, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END, min_mapq=0
    )
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
