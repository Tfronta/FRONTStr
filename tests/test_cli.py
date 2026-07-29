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


def test_call_command_is_gone() -> None:
    """LongTR is unwired: the command that ran it must not still be advertised."""
    result = CliRunner().invoke(app, ["--help"])
    assert "call" not in [line.split()[0] for line in result.stdout.splitlines() if line.strip()]
    assert CliRunner().invoke(app, ["call", "--help"]).exit_code != 0


def test_interpret_has_no_longtr_option() -> None:
    result = CliRunner().invoke(app, ["interpret", "--help"])
    assert result.exit_code == 0
    assert "longtr" not in result.stdout.lower()


def test_interpret_offers_the_log_flag() -> None:
    """The whole point of --log is discoverability; it must show in --help."""
    result = CliRunner().invoke(app, ["interpret", "--help"])
    assert result.exit_code == 0
    assert "--log" in result.stdout


def _option_names(command_name: str) -> set[str]:
    """Registered option strings for a subcommand.

    Introspected rather than scraped from ``--help``: Rich truncates the help
    to the terminal width, so a text assertion passes or fails depending on the
    window it runs in.
    """
    import typer.main

    group = typer.main.get_command(app)
    cmd = group.commands[command_name]  # type: ignore[attr-defined]
    return {opt for p in cmd.params for opt in getattr(p, "opts", [])}


def test_interpret_exposes_the_calling_thresholds() -> None:
    """ "Lower the min reads" has to be reachable without editing the panel."""
    opts = _option_names("interpret")
    for opt in ("--min-reads-third", "--min-phr", "--low-coverage-reads", "--balanced-ab-max"):
        assert opt in opts, f"{opt} missing from interpret"


def test_interpret_can_silence_the_parameter_echo() -> None:
    assert "--show-params" in _option_names("interpret")


def test_every_exposed_threshold_is_in_the_parameter_table() -> None:
    """A knob the CLI turns but the table cannot describe would be invisible in
    the report and unmarked when overridden."""
    from frontstr.params import PARAM_SPECS

    described = {s.name for s in PARAM_SPECS}
    exposed = {
        "--min-mapq": "min_mapq",
        "--flank-anchor": "flank_anchor",
        "--identity": "identity_threshold",
        "--len-tolerance": "len_tolerance_bp",
        "--analytical-thresh": "analytical_thresh",
        "--calling-thresh": "calling_thresh",
        "--min-phr": "min_phr_for_het",
        "--min-reads-third": "min_reads_third",
        "--low-coverage-reads": "low_coverage_reads",
        "--balanced-ab-max": "balanced_ab_max",
    }
    opts = _option_names("interpret")
    for flag, param in exposed.items():
        assert flag in opts, f"{flag} is no longer exposed"
        assert param in described, f"{param} is turned by {flag} but not in PARAM_SPECS"


def test_doctor_runs_without_a_bam() -> None:
    """The environment check has to be usable before you have data.

    The failure it exists for is quiet: with no POA backend FRONTStr still
    emits a full profile, from unpolished single reads, and the damage shows up
    as microvariants that are not in the sample.
    """
    result = CliRunner().invoke(app, ["doctor"])
    assert result.exit_code == 0, result.stdout
    assert "POA backend" in result.stdout
    assert "STRNaming" in result.stdout


def test_doctor_exposes_bed_and_region_options() -> None:
    assert "--bed" in _option_names("interpret")
    assert "--bed-coords" in _option_names("interpret")


def test_interpret_rejects_panel_and_bed_together(tmp_path: Path) -> None:
    bed = tmp_path / "r.bed"
    bed.write_text("chr11\t100\t200\tAATG\tTH01\n")
    result = CliRunner().invoke(
        app,
        ["interpret", "--bam", "x.bam", "--panel", "p.yaml", "--bed", str(bed)],
    )
    assert result.exit_code != 0


def test_interpret_requires_regions_from_somewhere() -> None:
    result = CliRunner().invoke(app, ["interpret", "--bam", "x.bam"])
    assert result.exit_code != 0
