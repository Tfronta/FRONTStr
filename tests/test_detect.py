"""Tests for input format detection."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from frontstr.errors import IngestError
from frontstr.ingest import detect_input
from frontstr.ingest.detect import InputKind


def test_detect_fastq_plain(tmp_fastq: Path) -> None:
    info = detect_input(tmp_fastq)
    assert info.kind == InputKind.FASTQ
    assert info.gzipped is False
    assert info.aligned is None
    assert info.size_bytes > 0


def test_detect_fastq_gz(tmp_fastq_gz: Path) -> None:
    info = detect_input(tmp_fastq_gz)
    assert info.kind == InputKind.FASTQ
    assert info.gzipped is True
    assert info.aligned is None


def test_detect_fasta_plain(tmp_fasta: Path) -> None:
    info = detect_input(tmp_fasta)
    assert info.kind == InputKind.FASTA
    assert info.gzipped is False


def test_detect_missing(tmp_path: Path) -> None:
    with pytest.raises(IngestError, match="does not exist"):
        detect_input(tmp_path / "nope.fastq")


def test_detect_empty(tmp_path: Path) -> None:
    empty = tmp_path / "empty.dat"
    empty.write_bytes(b"")
    with pytest.raises(IngestError, match="empty"):
        detect_input(empty)


def test_detect_garbage(tmp_path: Path) -> None:
    bad = tmp_path / "random.dat"
    bad.write_bytes(b"\x00\x01\x02\x03\x04\x05\x06\x07")
    with pytest.raises(IngestError, match="Unrecognised"):
        detect_input(bad)


def test_detect_gzipped_garbage(tmp_path: Path) -> None:
    bad = tmp_path / "weird.gz"
    with gzip.open(bad, "wb") as fh:
        fh.write(b"\x00\x01\x02\x03")
    with pytest.raises(IngestError, match="neither BAM"):
        detect_input(bad)


def test_detect_ignores_extension(tmp_path: Path) -> None:
    """Extension is irrelevant: a FASTQ named .bam must still be detected as FASTQ."""
    p = tmp_path / "mislabeled.bam"
    p.write_text("@r1\nACGT\n+\nIIII\n")
    info = detect_input(p)
    assert info.kind == InputKind.FASTQ


def test_detect_fake_bam_extension_gzip(tmp_path: Path) -> None:
    """A gzipped FASTQ with .bam extension is still FASTQ."""
    p = tmp_path / "fake.bam"
    with gzip.open(p, "wb") as fh:
        fh.write(b"@r1\nACGT\n+\nIIII\n")
    info = detect_input(p)
    assert info.kind == InputKind.FASTQ
    assert info.gzipped is True
