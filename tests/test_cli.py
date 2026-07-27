"""Tests for the FRONTStr CLI."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from frontstr.cli import app


def test_help() -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "Forensic" in result.stdout


def test_version() -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "FRONTStr" in result.stdout


def test_inspect_fastq(tmp_fastq: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["inspect", str(tmp_fastq)])
    assert result.exit_code == 0
    assert "fastq" in result.stdout.lower()


def test_inspect_missing() -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["inspect", "/nope/does/not/exist"])
    assert result.exit_code == 1


def test_run_not_implemented(tmp_path: Path) -> None:
    runner = CliRunner()
    p = tmp_path / "x.fastq"
    p.write_text("@r\nA\n+\nI\n")
    result = runner.invoke(
        app,
        [
            "run",
            "--input",
            str(p),
            "--sample",
            "S001",
            "--panel",
            str(tmp_path / "panel.yaml"),
            "--reference",
            str(tmp_path / "ref.fa"),
            "--out",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 64
