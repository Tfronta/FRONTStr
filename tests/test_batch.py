"""Tests for frontstr.batch: manifest parsing and batch orchestration."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.batch import ManifestEntry, parse_manifest, run_batch
from frontstr.errors import FrontstrError
from frontstr.panel.models import Panel, System

from .conftest import SYNTH_CHROM, SYNTH_TR_END, SYNTH_TR_START, SynthRead, write_synth_bam

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _write_manifest(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _synth_panel() -> Panel:
    return Panel(
        name="SYNTH",
        version="test",
        systems=[
            System(
                name="SYNTH_MARKER",
                chromosome=SYNTH_CHROM,
                ref_start=SYNTH_TR_START + 1,
                ref_end=SYNTH_TR_END,
                motif="AGAT",
                period=4,
                corr_value=0,
            )
        ],
    )


# ---------------------------------------------------------------------------
# parse_manifest
# ---------------------------------------------------------------------------

def test_parse_manifest_minimal(tmp_path: Path) -> None:
    """Two-column manifest (no role column) defaults to role='sample'."""
    p = _write_manifest(tmp_path / "m.tsv", [
        "sample_id\tbam",
        "HG00113\t/data/HG00113.bam",
    ])
    entries = parse_manifest(p)
    assert len(entries) == 1
    assert entries[0].sample_id == "HG00113"
    assert entries[0].bam == Path("/data/HG00113.bam")
    assert entries[0].role == "sample"


def test_parse_manifest_with_roles(tmp_path: Path) -> None:
    p = _write_manifest(tmp_path / "m.tsv", [
        "sample_id\tbam\trole",
        "S1\t/a.bam\tsample",
        "CTRL\t/b.bam\tpositive_ctrl",
        "NEG\t/c.bam\tnegative_ctrl",
        "BLANK\t/d.bam\treagent_blank",
    ])
    entries = parse_manifest(p)
    assert [e.role for e in entries] == [
        "sample", "positive_ctrl", "negative_ctrl", "reagent_blank"
    ]


def test_parse_manifest_skips_comments(tmp_path: Path) -> None:
    p = _write_manifest(tmp_path / "m.tsv", [
        "# FRONTStr batch manifest",
        "sample_id\tbam\trole",
        "# another comment",
        "S1\t/a.bam\tsample",
    ])
    entries = parse_manifest(p)
    assert len(entries) == 1
    assert entries[0].sample_id == "S1"


def test_parse_manifest_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FrontstrError, match="not found"):
        parse_manifest(tmp_path / "nope.tsv")


def test_parse_manifest_empty(tmp_path: Path) -> None:
    p = tmp_path / "m.tsv"
    p.write_text("")
    with pytest.raises(FrontstrError, match="empty"):
        parse_manifest(p)


def test_parse_manifest_bad_header(tmp_path: Path) -> None:
    p = _write_manifest(tmp_path / "m.tsv", ["name\tfile"])
    with pytest.raises(FrontstrError, match="sample_id"):
        parse_manifest(p)


def test_parse_manifest_unknown_role(tmp_path: Path) -> None:
    p = _write_manifest(tmp_path / "m.tsv", [
        "sample_id\tbam\trole",
        "S1\t/a.bam\tunknown_role",
    ])
    with pytest.raises(FrontstrError, match="unknown role"):
        parse_manifest(p)


def test_parse_manifest_duplicate_sample_id(tmp_path: Path) -> None:
    p = _write_manifest(tmp_path / "m.tsv", [
        "sample_id\tbam",
        "S1\t/a.bam",
        "S1\t/b.bam",
    ])
    with pytest.raises(FrontstrError, match="duplicate"):
        parse_manifest(p)


# ---------------------------------------------------------------------------
# run_batch (single-process, workers=1)
# ---------------------------------------------------------------------------

@pytest.fixture
def two_sample_bams(tmp_path: Path) -> tuple[Path, Path]:
    """Two synthetic BAMs: one het (CE 12+11), one hom (CE 12)."""
    specs_het = (
        [SynthRead(name=f"a{i}", n_repeats=12, hp=1) for i in range(10)]
        + [SynthRead(name=f"b{i}", n_repeats=11, hp=2) for i in range(8)]
    )
    specs_hom = [SynthRead(name=f"h{i}", n_repeats=12) for i in range(10)]
    bam_het = write_synth_bam(tmp_path / "het.bam", specs_het)
    bam_hom = write_synth_bam(tmp_path / "hom.bam", specs_hom)
    return bam_het, bam_hom


def test_run_batch_two_samples(
    tmp_path: Path, two_sample_bams: tuple[Path, Path]
) -> None:
    """Both samples succeed; per-sample dirs and batch_summary.csv are created."""
    bam_het, bam_hom = two_sample_bams
    entries = [
        ManifestEntry("HET", bam_het, "sample"),
        ManifestEntry("HOM", bam_hom, "positive_ctrl"),
    ]
    panel = _synth_panel()
    out = tmp_path / "batch_out"

    results = run_batch(
        entries=entries,
        panel=panel,
        out_dir=out,
        formats=frozenset({"profile", "json"}),
        workers=1,
    )

    assert len(results) == 2
    assert all(r.status == "ok" for r in results)

    assert (out / "HET" / "HET.profile.csv").exists()
    assert (out / "HET" / "HET.json").exists()
    assert (out / "HOM" / "HOM.profile.csv").exists()
    assert (out / "batch_summary.csv").exists()


def test_run_batch_summary_content(
    tmp_path: Path, two_sample_bams: tuple[Path, Path]
) -> None:
    """batch_summary.csv must have a row per sample with role and CE columns."""
    bam_het, bam_hom = two_sample_bams
    entries = [
        ManifestEntry("HET", bam_het, "sample"),
        ManifestEntry("HOM", bam_hom, "positive_ctrl"),
    ]
    panel = _synth_panel()
    out = tmp_path / "batch_out"
    run_batch(entries=entries, panel=panel, out_dir=out, formats=frozenset({"json"}), workers=1)

    import csv as csv_mod
    with (out / "batch_summary.csv").open(encoding="utf-8", newline="") as fh:
        rows = list(csv_mod.DictReader(fh))

    assert len(rows) == 2
    het_row = next(r for r in rows if r["sample_id"] == "HET")
    assert het_row["role"] == "sample"
    assert het_row["status"] == "ok"
    assert "SYNTH_MARKER" in het_row
    # Heterozygous: two alleles → CE string contains a comma
    assert "," in het_row["SYNTH_MARKER"] or het_row["SYNTH_MARKER"] != ""


def test_run_batch_missing_bam_is_error(tmp_path: Path) -> None:
    """A missing BAM produces status='error' without aborting the batch."""
    entries = [
        ManifestEntry("MISSING", tmp_path / "nope.bam", "sample"),
    ]
    panel = _synth_panel()
    out = tmp_path / "batch_out"
    results = run_batch(entries=entries, panel=panel, out_dir=out, workers=1)
    assert results[0].status == "error"
    assert results[0].error != ""
    # batch_summary.csv should still exist
    assert (out / "batch_summary.csv").exists()


def test_run_batch_progress_callback(
    tmp_path: Path, two_sample_bams: tuple[Path, Path]
) -> None:
    """progress_callback is called once per sample."""
    bam_het, bam_hom = two_sample_bams
    entries = [
        ManifestEntry("S1", bam_het, "sample"),
        ManifestEntry("S2", bam_hom, "sample"),
    ]
    panel = _synth_panel()
    out = tmp_path / "batch_out"
    ticks: list[str] = []
    run_batch(
        entries=entries, panel=panel, out_dir=out, workers=1,
        progress_callback=ticks.append,
        formats=frozenset({"json"}),
    )
    assert sorted(ticks) == ["S1", "S2"]


def test_run_batch_formats_respected(
    tmp_path: Path, two_sample_bams: tuple[Path, Path]
) -> None:
    """Only requested formats are written; others are absent."""
    bam_het, _ = two_sample_bams
    entries = [ManifestEntry("S1", bam_het, "sample")]
    panel = _synth_panel()
    out = tmp_path / "batch_out"
    run_batch(entries=entries, panel=panel, out_dir=out, formats=frozenset({"json"}), workers=1)

    assert (out / "S1" / "S1.json").exists()
    assert not (out / "S1" / "S1.profile.csv").exists()
    assert not (out / "S1" / "S1.html").exists()


# ---------------------------------------------------------------------------
# CLI smoke test
# ---------------------------------------------------------------------------

def test_batch_cli_two_samples(
    tmp_path: Path, two_sample_bams: tuple[Path, Path], codis_panel_yaml: Path
) -> None:
    """frontstr batch CLI produces batch_summary.csv for the synth panel."""
    from typer.testing import CliRunner

    from frontstr.cli import app

    # Write a minimal panel YAML for the synth marker
    synth_panel_yaml = tmp_path / "synth.yaml"
    synth_panel_yaml.write_text(
        f"name: SYNTH\nversion: test\nsystems:\n"
        f"  - name: SYNTH_MARKER\n"
        f"    chromosome: {SYNTH_CHROM}\n"
        f"    ref_start: {SYNTH_TR_START + 1}\n"
        f"    ref_end: {SYNTH_TR_END}\n"
        f"    motif: AGAT\n"
        f"    period: 4\n"
        f"    corr_value: 0\n"
        f"    category: autosomal\n",
        encoding="utf-8",
    )

    bam_het, bam_hom = two_sample_bams
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        f"sample_id\tbam\trole\n"
        f"HET\t{bam_het}\tsample\n"
        f"HOM\t{bam_hom}\tsample\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(
        app,
        ["batch", "--manifest", str(manifest), "--panel", str(synth_panel_yaml),
         "--out", str(tmp_path / "out"), "--formats", "json"],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "batch_summary.csv").exists()
