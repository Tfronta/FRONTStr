"""End-to-end demo of FRONTStr's evidence layer.

Builds a synthetic BAM with a known forensic profile (heterozygote CE 12/11
plus a couple of low-frequency reads at CE 10) and invokes the CLI to show
the per-allele integer coverage that FRONTStr produces.

Run:

    python examples/demo_evidence.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from tempfile import mkdtemp

THIS = Path(__file__).resolve()
sys.path.insert(0, str(THIS.parent.parent / "tests"))
from conftest import (  # noqa: E402  (import-after-sys-path-tweak intentional)
    SYNTH_CHROM,
    SynthRead,
    write_synth_bam,
)


def main() -> int:
    workdir = Path(mkdtemp(prefix="frontstr-demo-"))
    bam_path = workdir / "demo.bam"

    specs = (
        [SynthRead(name=f"a{i}", n_repeats=12, hp=1) for i in range(5)]
        + [SynthRead(name=f"b{i}", n_repeats=11, hp=2) for i in range(4)]
        + [SynthRead(name=f"c{i}", n_repeats=10) for i in range(2)]
    )
    write_synth_bam(bam_path, specs)
    print(f"synthetic BAM: {bam_path}")

    return subprocess.call(
        [
            sys.executable,
            "-m",
            "frontstr",
            "evidence",
            "--bam",
            str(bam_path),
            "--chrom",
            SYNTH_CHROM,
            "--start",
            "101",
            "--end",
            "148",
        ]
    )


if __name__ == "__main__":
    sys.exit(main())
