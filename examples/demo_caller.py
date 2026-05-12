"""Demo: parse a hand-crafted LongTR-style VCF via the `frontstr call` CLI.

Writes a tiny VCF that imitates LongTR's output (TH01 heterozygous + FGA with
a deletion + D3S1358 homozygous) and runs `frontstr call --parse-only` on it.

Run:  python examples/demo_caller.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from tempfile import mkdtemp

VCF_BODY = """##fileformat=VCFv4.2
##INFO=<ID=MOTIF,Number=.,Type=String,Description="Motif">
##INFO=<ID=PERIOD,Number=.,Type=String,Description="Period">
##INFO=<ID=BPDIFFS,Number=A,Type=Integer,Description="bp diff per ALT">
##INFO=<ID=INEXACT_ALLELE,Number=A,Type=Integer,Description="POA-derived ALT">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=Q,Number=1,Type=Float,Description="Posterior">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=PDP,Number=1,Type=String,Description="HP1|HP2 counts">
##FORMAT=<ID=ALLREADS,Number=1,Type=String,Description="bpdiff|reads">
##contig=<ID=chr3>
##contig=<ID=chr4>
##contig=<ID=chr11>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878
chr11\t2171100\tTH01\tAATG\tAATGAATG,AATGAATGAATG\t.\tPASS\tMOTIF=AATG;PERIOD=4;BPDIFFS=4,8;INEXACT_ALLELE=0,0\tGT:Q:DP:PDP:ALLREADS\t1/2:0.99:120:60|55:0|3;4|62;8|55
chr3\t45582240\tD3S1358\tTCTA\t.\t.\tPASS\tMOTIF=TCTA,TCTG;PERIOD=4,4\tGT:Q:DP:PDP:ALLREADS\t0/0:0.99:80:40|40:0|78
chr4\t154587740\tFGA\tAGAT\t<DEL>\t.\tPASS\tMOTIF=AGAT;PERIOD=4;BPDIFFS=-4;INEXACT_ALLELE=0\tGT:Q:DP:PDP:ALLREADS\t0/1:0.91:35:18|17:0|22;-4|13
"""


def main() -> int:
    workdir = Path(mkdtemp(prefix="frontstr-caller-demo-"))
    vcf = workdir / "demo.vcf"
    vcf.write_text(VCF_BODY)
    print(f"synthetic VCF: {vcf}")
    return subprocess.call(
        [
            sys.executable, "-m", "frontstr", "call",
            "--parse-only", str(vcf),
            "--bam", "/dev/null",
            "--panel", "/dev/null",
            "--reference", "/dev/null",
            "--out", str(workdir / "out"),
        ]
    )


if __name__ == "__main__":
    sys.exit(main())
