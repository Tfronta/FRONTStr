"""Tests for the LongTR runner.

Heavy lifting (subprocess) is skipped if the LongTR binary is not on PATH.
The argv builder is fully covered without the binary.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

from frontstr.caller.longtr import (
    LongTRInvocation,
    LongTRRunner,
    build_longtr_argv,
    run_longtr,
)
from frontstr.errors import CallerError
from frontstr.panel.models import Panel, System


@pytest.fixture
def fake_inputs(tmp_path: Path) -> dict[str, Path]:
    bam = tmp_path / "in.bam"
    bam.write_bytes(b"")
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 100 + "\n")
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t10\t40\tAGAT\tM1\n")
    return {"bam": bam, "fasta": fasta, "bed": bed, "vcf": tmp_path / "out.vcf.gz"}


def test_build_argv_ont_defaults(fake_inputs: dict[str, Path]) -> None:
    with patch.dict(os.environ, {"FRONTSTR_LONGTR_BIN": "echo"}):
        inv = build_longtr_argv(
            bam=fake_inputs["bam"],
            fasta=fake_inputs["fasta"],
            bed=fake_inputs["bed"],
            chrom=None,
            vcf_out=fake_inputs["vcf"],
            platform="ont",
        )
    assert inv.argv[0] == "echo"
    assert "--min-mean-qual" in inv.argv
    assert inv.argv[inv.argv.index("--min-mean-qual") + 1] == "10"
    assert "--alignment-params" in inv.argv
    assert "--phased-bam" not in inv.argv


def test_build_argv_hifi_qual(fake_inputs: dict[str, Path]) -> None:
    with patch.dict(os.environ, {"FRONTSTR_LONGTR_BIN": "echo"}):
        inv = build_longtr_argv(
            bam=fake_inputs["bam"],
            fasta=fake_inputs["fasta"],
            bed=fake_inputs["bed"],
            chrom=None,
            vcf_out=fake_inputs["vcf"],
            platform="hifi",
        )
    assert inv.argv[inv.argv.index("--min-mean-qual") + 1] == "30"
    assert "--alignment-params" not in inv.argv


def test_build_argv_extra_and_chrom(fake_inputs: dict[str, Path]) -> None:
    with patch.dict(os.environ, {"FRONTSTR_LONGTR_BIN": "echo"}):
        inv = build_longtr_argv(
            bam=fake_inputs["bam"],
            fasta=fake_inputs["fasta"],
            bed=fake_inputs["bed"],
            chrom="chr11",
            vcf_out=fake_inputs["vcf"],
            platform="ont",
            phased=True,
            skip_assembly=True,
            extra=["--use-unpaired"],
        )
    assert "--chrom" in inv.argv
    assert inv.argv[inv.argv.index("--chrom") + 1] == "chr11"
    assert "--phased-bam" in inv.argv
    assert "--skip-assembly" in inv.argv
    assert "--use-unpaired" in inv.argv


def test_build_argv_unknown_platform(fake_inputs: dict[str, Path]) -> None:
    with pytest.raises(CallerError, match="Unsupported platform"):
        build_longtr_argv(
            bam=fake_inputs["bam"],
            fasta=fake_inputs["fasta"],
            bed=fake_inputs["bed"],
            chrom=None,
            vcf_out=fake_inputs["vcf"],
            platform="illumina",
            binary="echo",
        )


def test_build_argv_missing_binary(
    fake_inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("FRONTSTR_LONGTR_BIN", raising=False)
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(CallerError, match="not found on PATH"):
        build_longtr_argv(
            bam=fake_inputs["bam"],
            fasta=fake_inputs["fasta"],
            bed=fake_inputs["bed"],
            chrom=None,
            vcf_out=fake_inputs["vcf"],
            platform="ont",
        )


def test_run_longtr_failure_propagates(tmp_path: Path) -> None:
    """A non-zero exit must raise CallerError with tail of the log."""
    inv = LongTRInvocation(
        argv=["false"],
        bam=tmp_path / "x.bam",
        bed=tmp_path / "x.bed",
        vcf=tmp_path / "x.vcf.gz",
    )
    with pytest.raises(CallerError, match="exited with code"):
        run_longtr(inv, log_path=tmp_path / "log.txt")


def test_run_longtr_missing_binary(tmp_path: Path) -> None:
    inv = LongTRInvocation(
        argv=["/this/binary/does/not/exist"],
        bam=tmp_path / "x.bam",
        bed=tmp_path / "x.bed",
        vcf=tmp_path / "x.vcf.gz",
    )
    with pytest.raises(CallerError, match="not found"):
        run_longtr(inv)


def test_run_longtr_success_but_no_vcf(tmp_path: Path) -> None:
    """A binary that exits 0 but produces no VCF is still an error."""
    inv = LongTRInvocation(
        argv=["true"],
        bam=tmp_path / "x.bam",
        bed=tmp_path / "x.bed",
        vcf=tmp_path / "missing.vcf.gz",
    )
    with pytest.raises(CallerError, match="did not write"):
        run_longtr(inv)


@pytest.mark.skipif(shutil.which("LongTR") is None, reason="LongTR not installed")
def test_runner_end_to_end(
    tmp_path: Path, fake_inputs: dict[str, Path]
) -> None:  # pragma: no cover
    """Integration test — only meaningful with a real LongTR install."""
    panel = Panel(
        name="t",
        version="0",
        systems=[
            System(name="M1", chromosome="chr1", ref_start=10, ref_end=40, motif="AGAT", period=4)
        ],
    )
    runner = LongTRRunner(panel=panel, reference=fake_inputs["fasta"], platform="ont")
    # Just exercise argv construction; we expect an error because the BAM is empty.
    with pytest.raises(CallerError):
        runner.run(bam=fake_inputs["bam"], out_dir=tmp_path / "out")
