"""Shared pytest fixtures.

The synthetic BAM fixtures build minimal in-memory BAMs that exercise
:mod:`frontstr.evidence` without requiring real ONT data. The reference is
``chr1`` with the TR at 0-based ``[100, 148)`` (1-based inclusive 101-148),
motif ``AGAT``, ``corr_value=0``, so:

  - 12 repeats → 48 bp TR, CE = 12.0  (the "reference" allele)
  - 11 repeats → 44 bp TR, CE = 11.0
  - 10 repeats → 40 bp TR, CE = 10.0
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pysam
import pytest


@pytest.fixture
def tmp_fastq(tmp_path: Path) -> Path:
    p = tmp_path / "sample.fastq"
    p.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGTACGT\n+\nIIIIIIII\n")
    return p


@pytest.fixture
def tmp_fastq_gz(tmp_path: Path) -> Path:
    p = tmp_path / "sample.fastq.gz"
    with gzip.open(p, "wb") as fh:
        fh.write(b"@r1\nACGT\n+\nIIII\n")
    return p


@pytest.fixture
def tmp_fasta(tmp_path: Path) -> Path:
    p = tmp_path / "sample.fasta"
    p.write_text(">r1\nACGT\n>r2\nACGTACGT\n")
    return p


@pytest.fixture
def codis_panel_yaml() -> Path:
    """Path to the example CODIS 20 panel shipped with the repo."""
    here = Path(__file__).resolve().parent
    return here.parent / "examples" / "panels" / "codis_20_grch38.yaml"


@pytest.fixture
def demo_catalog_json() -> Path:
    """Path to the demonstration allele catalog shipped with the repo."""
    here = Path(__file__).resolve().parent
    return here.parent / "examples" / "catalogs" / "demo_seed.json"


# ---------------------------------------------------------------------------
# Synthetic BAM helpers for evidence-layer tests
# ---------------------------------------------------------------------------

# Locus coordinates (0-based half-open)
SYNTH_CHROM = "chr1"
SYNTH_TR_START = 100
SYNTH_TR_END = 148  # 12 repeats of "AGAT" = 48 bp
SYNTH_MOTIF = "AGAT"
SYNTH_FLANK_LEN = 100
SYNTH_CHR_LEN = 1000


@dataclass(slots=True)
class SynthRead:
    """Specification for one synthetic read in the test BAM."""

    name: str
    n_repeats: int
    hp: int | None = None
    mapq: int = 60
    reverse: bool = False
    substitution_at: int | None = None  # position inside the TR for a single SNV


def _build_read(spec: SynthRead, header: pysam.AlignmentHeader) -> pysam.AlignedSegment:
    """Create one aligned read with deterministic CIGAR encoding length changes."""
    flank_l = "A" * SYNTH_FLANK_LEN
    flank_r = "T" * SYNTH_FLANK_LEN
    tr = SYNTH_MOTIF * spec.n_repeats
    if spec.substitution_at is not None and 0 <= spec.substitution_at < len(tr):
        original = tr[spec.substitution_at]
        replacement = "G" if original != "G" else "C"
        tr = tr[: spec.substitution_at] + replacement + tr[spec.substitution_at + 1 :]
    seq = flank_l + tr + flank_r
    quals = "I" * len(seq)

    a = pysam.AlignedSegment(header=header)
    a.query_name = spec.name
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array(quals)
    a.reference_id = 0
    a.reference_start = 0
    a.mapping_quality = spec.mapq
    a.flag = 16 if spec.reverse else 0  # reverse strand bit, no other flags

    ref_tr_len = (SYNTH_TR_END - SYNTH_TR_START)  # 48
    read_tr_len = 4 * spec.n_repeats
    diff = ref_tr_len - read_tr_len
    if diff == 0:
        a.cigartuples = [(0, SYNTH_FLANK_LEN), (0, read_tr_len), (0, SYNTH_FLANK_LEN)]
    elif diff > 0:
        # Read is shorter → represent missing bp as deletion at end of TR
        a.cigartuples = [
            (0, SYNTH_FLANK_LEN),
            (0, read_tr_len),
            (2, diff),
            (0, SYNTH_FLANK_LEN),
        ]
    else:
        # Read is longer → insertion of |diff| bp at end of TR
        a.cigartuples = [
            (0, SYNTH_FLANK_LEN),
            (0, ref_tr_len),
            (1, -diff),
            (0, SYNTH_FLANK_LEN),
        ]
    a.set_tag("RG", "rg1")
    if spec.hp is not None:
        a.set_tag("HP", spec.hp)
    return a


def write_synth_bam(out_path: Path, specs: Iterable[SynthRead]) -> Path:
    """Write a synthetic, indexed, coordinate-sorted BAM at ``out_path``.

    Returns the path to the indexed BAM.
    """
    header_dict = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": SYNTH_CHROM, "LN": SYNTH_CHR_LEN}],
        "RG": [{"ID": "rg1", "SM": "synthetic", "LB": "lib1", "PL": "ONT"}],
    }
    raw = out_path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(raw), "wb", header=header_dict) as bam:
        for spec in specs:
            bam.write(_build_read(spec, bam.header))
    pysam.sort("-o", str(out_path), str(raw))
    pysam.index(str(out_path))
    raw.unlink(missing_ok=True)
    return out_path


@pytest.fixture
def synth_bam_heterozygous(tmp_path: Path) -> Path:
    """11 reads: 5 x CE12 (HP1), 4 x CE11 (HP2), 2 x CE10 (no HP)."""
    specs: list[SynthRead] = []
    for i in range(5):
        specs.append(SynthRead(name=f"a{i}", n_repeats=12, hp=1))
    for i in range(4):
        specs.append(SynthRead(name=f"b{i}", n_repeats=11, hp=2))
    for i in range(2):
        specs.append(SynthRead(name=f"c{i}", n_repeats=10, hp=None))
    return write_synth_bam(tmp_path / "het.bam", specs)


@pytest.fixture
def synth_bam_homozygous(tmp_path: Path) -> Path:
    """8 reads all CE12, mixed strands, no HP tags."""
    specs = [
        SynthRead(name=f"r{i}", n_repeats=12, reverse=(i % 2 == 0))
        for i in range(8)
    ]
    return write_synth_bam(tmp_path / "hom.bam", specs)


@pytest.fixture
def synth_bam_with_noise(tmp_path: Path) -> Path:
    """5 reads CE12 + 1 read CE12 with one substitution.

    Designed to test identity-threshold behaviour:
    - identity=0.97 → 1 cluster of 6
    - identity=0.99 → 2 clusters (5 + 1)
    """
    specs: list[SynthRead] = []
    for i in range(5):
        specs.append(SynthRead(name=f"clean{i}", n_repeats=12))
    specs.append(SynthRead(name="noise", n_repeats=12, substitution_at=20))
    return write_synth_bam(tmp_path / "noise.bam", specs)


@pytest.fixture
def synth_bam_low_mapq(tmp_path: Path) -> Path:
    """3 high-MAPQ reads + 3 low-MAPQ reads at the same locus."""
    specs = [
        SynthRead(name="hi0", n_repeats=12, mapq=60),
        SynthRead(name="hi1", n_repeats=12, mapq=60),
        SynthRead(name="hi2", n_repeats=12, mapq=60),
        SynthRead(name="lo0", n_repeats=10, mapq=5),
        SynthRead(name="lo1", n_repeats=10, mapq=5),
        SynthRead(name="lo2", n_repeats=10, mapq=5),
    ]
    return write_synth_bam(tmp_path / "mapq.bam", specs)
