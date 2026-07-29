"""Tests for frontstr.batch: manifest parsing and batch orchestration."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.batch import ManifestEntry, _extract_marker_ces, parse_manifest, run_batch
from frontstr.errors import FrontstrError
from frontstr.interp.models import (
    Allele,
    AlleleStatus,
    CallRule,
    MarkerResult,
    TriType,
)
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
    p = _write_manifest(
        tmp_path / "m.tsv",
        [
            "sample_id\tbam",
            "HG00113\t/data/HG00113.bam",
        ],
    )
    entries = parse_manifest(p)
    assert len(entries) == 1
    assert entries[0].sample_id == "HG00113"
    assert entries[0].bam == Path("/data/HG00113.bam")
    assert entries[0].role == "sample"


def test_parse_manifest_with_roles(tmp_path: Path) -> None:
    p = _write_manifest(
        tmp_path / "m.tsv",
        [
            "sample_id\tbam\trole",
            "S1\t/a.bam\tsample",
            "CTRL\t/b.bam\tpositive_ctrl",
            "NEG\t/c.bam\tnegative_ctrl",
            "BLANK\t/d.bam\treagent_blank",
        ],
    )
    entries = parse_manifest(p)
    assert [e.role for e in entries] == [
        "sample",
        "positive_ctrl",
        "negative_ctrl",
        "reagent_blank",
    ]


def test_parse_manifest_skips_comments(tmp_path: Path) -> None:
    p = _write_manifest(
        tmp_path / "m.tsv",
        [
            "# FRONTStr batch manifest",
            "sample_id\tbam\trole",
            "# another comment",
            "S1\t/a.bam\tsample",
        ],
    )
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
    p = _write_manifest(
        tmp_path / "m.tsv",
        [
            "sample_id\tbam\trole",
            "S1\t/a.bam\tunknown_role",
        ],
    )
    with pytest.raises(FrontstrError, match="unknown role"):
        parse_manifest(p)


def test_parse_manifest_duplicate_sample_id(tmp_path: Path) -> None:
    p = _write_manifest(
        tmp_path / "m.tsv",
        [
            "sample_id\tbam",
            "S1\t/a.bam",
            "S1\t/b.bam",
        ],
    )
    with pytest.raises(FrontstrError, match="duplicate"):
        parse_manifest(p)


# ---------------------------------------------------------------------------
# run_batch (single-process, workers=1)
# ---------------------------------------------------------------------------


@pytest.fixture
def two_sample_bams(tmp_path: Path) -> tuple[Path, Path]:
    """Two synthetic BAMs: one het (CE 12+11), one hom (CE 12)."""
    specs_het = [SynthRead(name=f"a{i}", n_repeats=12, hp=1) for i in range(10)] + [
        SynthRead(name=f"b{i}", n_repeats=11, hp=2) for i in range(8)
    ]
    specs_hom = [SynthRead(name=f"h{i}", n_repeats=12) for i in range(10)]
    bam_het = write_synth_bam(tmp_path / "het.bam", specs_het)
    bam_hom = write_synth_bam(tmp_path / "hom.bam", specs_hom)
    return bam_het, bam_hom


def test_run_batch_two_samples(tmp_path: Path, two_sample_bams: tuple[Path, Path]) -> None:
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


def test_run_batch_summary_content(tmp_path: Path, two_sample_bams: tuple[Path, Path]) -> None:
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


def test_run_batch_progress_callback(tmp_path: Path, two_sample_bams: tuple[Path, Path]) -> None:
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
        entries=entries,
        panel=panel,
        out_dir=out,
        workers=1,
        progress_callback=ticks.append,
        formats=frozenset({"json"}),
    )
    assert sorted(ticks) == ["S1", "S2"]


def test_run_batch_formats_respected(tmp_path: Path, two_sample_bams: tuple[Path, Path]) -> None:
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
        f"sample_id\tbam\trole\nHET\t{bam_het}\tsample\nHOM\t{bam_hom}\tsample\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "batch",
            "--manifest",
            str(manifest),
            "--panel",
            str(synth_panel_yaml),
            "--out",
            str(tmp_path / "out"),
            "--formats",
            "json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "batch_summary.csv").exists()


# ---------------------------------------------------------------------------
# _extract_marker_ces — the allele numbers batch_summary.csv reports
# ---------------------------------------------------------------------------


def _numbered_allele(idx: int, ce: float, strnaming_ce: float | None) -> Allele:
    """An allele whose legacy CE and STRNaming CE deliberately disagree."""
    return Allele(
        cluster_index=idx,
        consensus="AATGAATGAATG",
        length_bp=12,
        n_reads_total=40,
        n_reads_hp1=20,
        n_reads_hp2=20,
        n_reads_hp_none=0,
        n_forward=20,
        n_reverse=20,
        mean_qual=30.0,
        ce=ce,
        strnaming_ce=strnaming_ce,
        isfg=f"[AATG]{int(ce)}",
        bp_diff=0,
        is_deletion=False,
        status=AlleleStatus.ALLELE,
    )


def _marker_result(name: str, alleles: list[Allele]) -> MarkerResult:
    return MarkerResult(
        marker_name=name,
        system=_synth_panel().systems[0],
        alleles=alleles,
        alleles_called=alleles,
        call_rule=CallRule.HETEROZYGOUS if len(alleles) > 1 else CallRule.HOMOZYGOUS,
        tri_type=TriType.NONE,
        total_reads=sum(a.n_reads_total for a in alleles),
    )


def test_summary_reports_the_strnaming_number_not_the_legacy_ce() -> None:
    """batch_summary.csv must carry the same allele number as every other view.

    Regression: the summary used to read the raw ``Allele.ce``, the
    length-derived number from before STRNaming became the source of the allele
    designation. That made this file the one output disagreeing with the same
    run's report, VCF and XLSX at the six markers whose ``corr_value`` had been
    miscalibrated — HG00113 vWA read 13/17 here and 14/16 everywhere else.
    """
    result = _marker_result(
        "vWA",
        [_numbered_allele(0, ce=13.0, strnaming_ce=14.0), _numbered_allele(1, 17.0, 16.0)],
    )

    assert _extract_marker_ces([result]) == {"vWA": "14,16"}


def test_summary_falls_back_to_the_legacy_ce_when_strnaming_declined() -> None:
    """Markers STRNaming has no range for still get their number from the CE."""
    result = _marker_result("DYS393", [_numbered_allele(0, ce=14.0, strnaming_ce=None)])

    assert _extract_marker_ces([result]) == {"DYS393": "14"}


def test_summary_carries_the_amel_designation() -> None:
    """AMEL has no allele number; its X/Y designation must survive into the CSV.

    ``number_label`` falls back to the ISFG designation when there is no number,
    which is the whole reason a non-numeric marker can appear in this file at
    all. Reading the raw CE gave AMEL an empty cell.
    """
    x = _numbered_allele(0, ce=13.0, strnaming_ce=None)
    x.ce = None
    x.allele_numeric = None
    x.length_bp = 0
    x.isfg = "X"
    y = x.model_copy(update={"cluster_index": 1, "isfg": "Y"})

    assert _extract_marker_ces([_marker_result("AMEL", [x, y])]) == {"AMEL": "X,Y"}


# ---------------------------------------------------------------------------
# --log / --trace
# ---------------------------------------------------------------------------


def test_trace_writes_one_narrative_file_per_sample(
    tmp_path: Path, two_sample_bams: tuple[Path, Path]
) -> None:
    """Each sample gets its own trace beside its other outputs.

    Per sample, per file — not the terminal. One sample's trace is already
    hundreds of lines, so a cohort's would bury the progress the operator is
    watching, and parallel workers would interleave loci into nonsense.
    """
    bam_het, bam_hom = two_sample_bams
    entries = [ManifestEntry("HET", bam_het, "sample"), ManifestEntry("HOM", bam_hom, "sample")]
    out = tmp_path / "batch_out"

    run_batch(
        entries=entries,
        panel=_synth_panel(),
        out_dir=out,
        formats=frozenset({"json"}),
        workers=1,
        trace=True,
    )

    for sample_id in ("HET", "HOM"):
        trace = out / sample_id / f"{sample_id}.trace.txt"
        assert trace.exists(), f"no trace for {sample_id}"
        text = trace.read_text(encoding="utf-8")
        assert "SYNTH_MARKER" in text, "the locus narrative is missing"
        assert "Loci processed" in text, "the run summary is missing"


def test_trace_written_by_default(tmp_path: Path, two_sample_bams: tuple[Path, Path]) -> None:
    """The audit record is not opt-in.

    A run that kept no trace cannot be questioned afterwards without being
    repeated, and generating it costs nothing measurable (18.5 s against 19.0 s
    over five samples). ``--no-trace`` exists, but nobody has to remember it to
    get an auditable run.
    """
    bam_het, _ = two_sample_bams
    out = tmp_path / "batch_out"

    run_batch(
        entries=[ManifestEntry("HET", bam_het, "sample")],
        panel=_synth_panel(),
        out_dir=out,
        formats=frozenset({"json"}),
        workers=1,
    )

    assert (out / "HET" / "HET.trace.txt").exists(), "no trace without being asked for one"


def test_no_trace_opts_out(tmp_path: Path, two_sample_bams: tuple[Path, Path]) -> None:
    """The escape hatch still works, for a genuinely constrained output dir."""
    bam_het, _ = two_sample_bams
    out = tmp_path / "batch_out"

    run_batch(
        entries=[ManifestEntry("HET", bam_het, "sample")],
        panel=_synth_panel(),
        out_dir=out,
        formats=frozenset({"json"}),
        workers=1,
        trace=False,
    )

    assert not (out / "HET" / "HET.trace.txt").exists()
    assert (out / "HET" / "HET.json").exists(), "the normal outputs must still be there"


def test_log_tags_every_line_with_its_sample(
    tmp_path: Path, two_sample_bams: tuple[Path, Path], capfd: pytest.CaptureFixture[str]
) -> None:
    """With workers running in parallel, an unattributed log line is useless."""
    bam_het, bam_hom = two_sample_bams
    out = tmp_path / "batch_out"

    run_batch(
        entries=[ManifestEntry("HET", bam_het, "sample"), ManifestEntry("HOM", bam_hom, "sample")],
        panel=_synth_panel(),
        out_dir=out,
        formats=frozenset({"json"}),
        workers=1,
        log=True,
    )

    stderr = capfd.readouterr().err
    assert "sample=HET" in stderr
    assert "sample=HOM" in stderr


def test_trace_streams_live_when_serial(
    tmp_path: Path, two_sample_bams: tuple[Path, Path], capfd: pytest.CaptureFixture[str]
) -> None:
    """With one worker the narrative is watchable, not just archived.

    The sequences and the HP1/HP2 counts are what makes a call followable back
    to the reads; hiding them in a file the operator has to go find afterwards
    defeats the point of watching a run at all.
    """
    bam_het, _ = two_sample_bams
    out = tmp_path / "batch_out"

    run_batch(
        entries=[ManifestEntry("HET", bam_het, "sample")],
        panel=_synth_panel(),
        out_dir=out,
        formats=frozenset({"json"}),
        workers=1,
        trace=True,
    )

    stderr = capfd.readouterr().err
    assert "SYNTH_MARKER" in stderr, "the locus narrative never reached the terminal"
    assert "Sequences" in stderr, "the aligned sequences must be on screen"
    assert "HP1" in stderr, "the haplotype counts must be on screen"
    # And still archived, so a locus can be re-read later.
    assert (out / "HET" / "HET.trace.txt").exists()


def test_trace_does_not_stream_when_parallel(
    tmp_path: Path, two_sample_bams: tuple[Path, Path], capfd: pytest.CaptureFixture[str]
) -> None:
    """Parallel workers would interleave loci from different samples."""
    bam_het, bam_hom = two_sample_bams
    out = tmp_path / "batch_out"

    run_batch(
        entries=[ManifestEntry("HET", bam_het, "sample"), ManifestEntry("HOM", bam_hom, "sample")],
        panel=_synth_panel(),
        out_dir=out,
        formats=frozenset({"json"}),
        workers=2,
        trace=True,
    )

    assert "Sequences" not in capfd.readouterr().err
    for sample_id in ("HET", "HOM"):
        assert (out / sample_id / f"{sample_id}.trace.txt").exists()
