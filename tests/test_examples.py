"""Smoke test: every shipped example must actually run.

``examples/`` is in the sdist and is the first thing a new user runs, but
nothing executed it. Four of the five demos had been passing ``--longtr-vcf``
to the CLI for as long as that option had not existed, so each one exited 2 on
a usage error. Lint cannot catch that: the flag is a string in an argv list.

Each demo builds a synthetic BAM in a temp directory and shells out to the CLI,
so this is an end-to-end check of the command line as a user meets it. The
whole file runs in about five seconds.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

DEMOS = sorted(p.name for p in EXAMPLES.glob("demo_*.py"))


def test_the_demos_are_discovered() -> None:
    """Guard against the glob silently matching nothing after a rename."""
    assert DEMOS, f"no demo_*.py found under {EXAMPLES}"


@pytest.mark.parametrize("demo", DEMOS)
def test_demo_runs(demo: str) -> None:
    proc = subprocess.run(
        [sys.executable, str(EXAMPLES / demo)],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, f"{demo} exited {proc.returncode}\n{proc.stdout}\n{proc.stderr}"
