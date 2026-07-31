"""Demo: export all formats from a synthetic FRONTStr run.

Builds the same demo BAM and panel as ``demo_interpret.py`` then calls::

    python -m frontstr export --formats profile,evidence,seqs,json,html

…and lists everything that was written.

Run:  python examples/demo_export.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from tempfile import mkdtemp

THIS = Path(__file__).resolve()
sys.path.insert(0, str(THIS.parent))
from demo_interpret import (  # noqa: E402
    PANEL_YAML,
    _build_bam,
)


def main() -> int:
    workdir = Path(mkdtemp(prefix="frontstr-export-demo-"))
    print(f"workdir: {workdir}")
    panel = workdir / "panel.yaml"
    panel.write_text(PANEL_YAML)
    bam = _build_bam(workdir / "demo.bam")
    out_dir = workdir / "exports"
    return subprocess.call(
        [
            sys.executable, "-m", "frontstr", "export",
            "--bam", str(bam),
            "--panel", str(panel),
            "--sample", "FRONTStr-Demo",
            "--operator", "demo@frontstr",
            "--out-dir", str(out_dir),
            "--formats", "profile,evidence,seqs,json,html",
        ]
    )


if __name__ == "__main__":
    sys.exit(main())
