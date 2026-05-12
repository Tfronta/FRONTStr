"""Demo: generate a full FRONTStr HTML report on a synthetic dataset.

Reuses the ``demo_interpret.py`` synthetic-BAM builder, then invokes::

    python -m frontstr report --bam ... --panel ... --longtr-vcf ... --out report.html

…and prints the path so you can open it in a browser.

Run:  python examples/demo_report.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from tempfile import mkdtemp

THIS = Path(__file__).resolve()
sys.path.insert(0, str(THIS.parent))
from demo_interpret import (  # noqa: E402
    LONGTR_VCF,
    PANEL_YAML,
    _build_bam,
)


def main() -> int:
    workdir = Path(mkdtemp(prefix="frontstr-report-demo-"))
    print(f"workdir: {workdir}")
    panel = workdir / "panel.yaml"
    panel.write_text(PANEL_YAML)
    vcf = workdir / "longtr.vcf"
    vcf.write_text(LONGTR_VCF)
    bam = _build_bam(workdir / "demo.bam")
    out_html = workdir / "report.html"
    rc = subprocess.call(
        [
            sys.executable, "-m", "frontstr", "report",
            "--bam", str(bam),
            "--panel", str(panel),
            "--longtr-vcf", str(vcf),
            "--sample", "FRONTStr-Demo",
            "--operator", "demo@frontstr",
            "--run-id", "demo-001",
            "--out", str(out_html),
        ]
    )
    if rc == 0:
        size_kb = out_html.stat().st_size / 1024
        print(f"\nHTML report ready: file://{out_html}  ({size_kb:.1f} KB)")
    return rc


if __name__ == "__main__":
    sys.exit(main())
