"""Demo: full FRONTStr pipeline end-to-end on a synthetic dataset.

We build:
- A synthetic BAM with TWO loci stitched in: a heterozygous TH01-like locus
  and a TPOX-like triallelic Type II locus on the same chr1.
- A minimal panel YAML covering both loci.
- A LongTR-style VCF asserting one of the calls.

Then we run::

    python -m frontstr interpret --bam ... --panel ... --longtr-vcf ...

…and print the table of forensic calls.

Run:  python examples/demo_interpret.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from tempfile import mkdtemp

import pysam

PANEL_YAML = """name: demo-panel
version: "0.1"
reference_build: GRCh38
description: Demo two-locus panel for FRONTStr Phase 1.5
systems:
  - name: TH01_LIKE
    chromosome: chr1
    ref_start: 101
    ref_end: 148
    motif: AGAT
    period: 4
    corr_value: 0
    category: autosomal
  - name: TPOX_LIKE
    chromosome: chr1
    ref_start: 301
    ref_end: 332
    motif: AATG
    period: 4
    corr_value: 0
    category: autosomal
    allow_triallelic: true
    tri_balanced_thr: 0.5
"""

LONGTR_VCF = """##fileformat=VCFv4.2
##INFO=<ID=MOTIF,Number=.,Type=String,Description="Motif">
##INFO=<ID=PERIOD,Number=.,Type=String,Description="Period">
##INFO=<ID=BPDIFFS,Number=A,Type=Integer,Description="bp diff per ALT">
##INFO=<ID=INEXACT_ALLELE,Number=A,Type=Integer,Description="Inexact ALT">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=Q,Number=1,Type=Float,Description="Q">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=PDP,Number=1,Type=String,Description="PDP">
##FORMAT=<ID=ALLREADS,Number=1,Type=String,Description="ALLREADS">
##contig=<ID=chr1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdemo
chr1\t101\tTH01_LIKE\tAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT\tAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT\t.\tPASS\tMOTIF=AGAT;PERIOD=4;BPDIFFS=-4;INEXACT_ALLELE=0\tGT:Q:DP:PDP:ALLREADS\t0/1:0.99:58:30|25:0|33;-4|25
"""


def _build_chr1_fasta(out_path: Path) -> None:
    """Build a 1 kb chr1 with two TR regions hand-placed at the right offsets."""
    flank_l = "A" * 100
    th01_ref = "AGAT" * 12          # positions 100..148 (0-based), 1-based 101..148
    middle = "T" * (300 - len(flank_l) - len(th01_ref))  # filler
    tpox_ref = "AATG" * 8           # 32 bp, 1-based 301..332
    tail = "G" * (1000 - len(flank_l) - len(th01_ref) - len(middle) - len(tpox_ref))
    seq = flank_l + th01_ref + middle + tpox_ref + tail
    assert len(seq) == 1000
    out_path.write_text(f">chr1\n{seq}\n")


def _add_read(bam: pysam.AlignmentFile, *, name: str, start: int, seq: str,
              ref_len: int, mapq: int = 60, hp: int | None = None,
              reverse: bool = False) -> None:
    a = pysam.AlignedSegment(header=bam.header)
    a.query_name = name
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = mapq
    a.flag = 16 if reverse else 0
    # Construct CIGAR: flank-match, tr-match, deletion or insertion, flank-match
    flank_l = 100
    flank_r = len(seq) - flank_l - (len(seq) - flank_l - flank_l)  # placeholder
    # Simpler: assume seq = flank_l*A + TR + flank_l*T; len(TR)=read_tr_len, ref_tr_len=ref_len
    read_tr_len = len(seq) - 2 * flank_l
    diff = ref_len - read_tr_len
    if diff == 0:
        a.cigartuples = [(0, flank_l), (0, read_tr_len), (0, flank_l)]
    elif diff > 0:
        a.cigartuples = [(0, flank_l), (0, read_tr_len), (2, diff), (0, flank_l)]
    else:
        a.cigartuples = [(0, flank_l), (0, ref_len), (1, -diff), (0, flank_l)]
    a.set_tag("RG", "rg1")
    if hp is not None:
        a.set_tag("HP", hp)
    bam.write(a)


def _make_read(motif: str, n_repeats: int, ref_len: int, flank_l_seq: str,
               flank_r_seq: str) -> str:
    """flank_l_seq and flank_r_seq are pre-computed flanking sequences."""
    return flank_l_seq + motif * n_repeats + flank_r_seq


def _build_bam(out_bam: Path) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "RG": [{"ID": "rg1", "SM": "demo", "LB": "lib1", "PL": "ONT"}],
    }
    raw = out_bam.with_suffix(".unsorted.bam")
    flank_left_th01 = "A" * 100
    flank_right_th01 = "T" * 100
    flank_left_tpox = "C" * 100
    flank_right_tpox = "G" * 100

    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        # TH01_LIKE: 30x CE12 (HP1) + 25x CE11 (HP2) + 3x CE10 (stutter)
        for i in range(30):
            seq = _make_read("AGAT", 12, 48, flank_left_th01, flank_right_th01)
            _add_read(bam, name=f"th01_a{i}", start=0, seq=seq, ref_len=48, hp=1)
        for i in range(25):
            seq = _make_read("AGAT", 11, 48, flank_left_th01, flank_right_th01)
            _add_read(bam, name=f"th01_b{i}", start=0, seq=seq, ref_len=48, hp=2)
        for i in range(3):
            seq = _make_read("AGAT", 10, 48, flank_left_th01, flank_right_th01)
            _add_read(bam, name=f"th01_c{i}", start=0, seq=seq, ref_len=48)
        # TPOX_LIKE: balanced triallelic Type II - 22x CE8 + 20x CE9 + 20x CE7
        for i in range(22):
            seq = _make_read("AATG", 8, 32, flank_left_tpox, flank_right_tpox)
            _add_read(bam, name=f"tpox_a{i}", start=200, seq=seq, ref_len=32, hp=1)
        for i in range(20):
            seq = _make_read("AATG", 9, 32, flank_left_tpox, flank_right_tpox)
            _add_read(bam, name=f"tpox_b{i}", start=200, seq=seq, ref_len=32, hp=2)
        for i in range(20):
            seq = _make_read("AATG", 7, 32, flank_left_tpox, flank_right_tpox)
            _add_read(bam, name=f"tpox_c{i}", start=200, seq=seq, ref_len=32)
    pysam.sort("-o", str(out_bam), str(raw))
    pysam.index(str(out_bam))
    raw.unlink(missing_ok=True)
    return out_bam


def main() -> int:
    workdir = Path(mkdtemp(prefix="frontstr-interp-demo-"))
    print(f"workdir: {workdir}")
    panel = workdir / "panel.yaml"
    panel.write_text(PANEL_YAML)
    vcf = workdir / "longtr.vcf"
    vcf.write_text(LONGTR_VCF)
    bam = _build_bam(workdir / "demo.bam")
    print(f"synthetic BAM: {bam}")
    return subprocess.call(
        [
            sys.executable, "-m", "frontstr", "interpret",
            "--bam", str(bam),
            "--panel", str(panel),
            "--longtr-vcf", str(vcf),
        ]
    )


if __name__ == "__main__":
    sys.exit(main())
