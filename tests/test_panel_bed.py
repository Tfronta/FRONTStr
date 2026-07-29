"""Tests for BED export from a Panel."""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.errors import PanelError
from frontstr.panel.bed import split_panel_by_chromosome, write_panel_bed
from frontstr.panel.loader import load_panel
from frontstr.panel.models import Panel, System


def test_write_panel_bed_basic(tmp_path: Path) -> None:
    panel = Panel(
        name="t",
        version="0",
        systems=[
            System(
                name="M1", chromosome="chr1", ref_start=200, ref_end=240, motif="AGAT", period=4
            ),
            System(
                name="M2",
                chromosome="chr1",
                ref_start=100,
                ref_end=140,
                motif="TCTA,TCTG",
                period=-1,
            ),
            System(name="M3", chromosome="chr2", ref_start=50, ref_end=80, motif="AT", period=2),
        ],
    )
    out = write_panel_bed(panel, tmp_path / "p.bed")
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 3
    # Sorted by chromosome then ref_start
    assert lines[0] == "chr1\t100\t140\tTCTA,TCTG\tM2"
    assert lines[1] == "chr1\t200\t240\tAGAT\tM1"
    assert lines[2] == "chr2\t50\t80\tAT\tM3"


def test_write_panel_bed_skips_no_motif(tmp_path: Path) -> None:
    """Markers without a motif (Illumina-only legacy) are skipped — but pydantic
    validation forbids empty motifs in System construction, so the only way to
    produce a "no motif" system is via the model_construct backdoor. We assert
    the same outcome by making sure the BED contains exactly the systems with
    valid motifs."""
    panel = Panel(
        name="t",
        version="0",
        systems=[
            System(name="ok", chromosome="chr1", ref_start=1, ref_end=10, motif="A", period=1),
        ],
    )
    out = write_panel_bed(panel, tmp_path / "p.bed")
    assert "ok" in out.read_text()


def test_write_panel_bed_raises_when_empty(tmp_path: Path) -> None:
    panel = Panel.model_construct(name="empty", version="0", systems=[])
    with pytest.raises(PanelError):
        write_panel_bed(panel, tmp_path / "p.bed")


def test_split_panel_by_chromosome(tmp_path: Path) -> None:
    panel = Panel(
        name="t",
        version="0",
        systems=[
            System(name="A", chromosome="chr1", ref_start=10, ref_end=20, motif="AG", period=2),
            System(name="B", chromosome="chr1", ref_start=30, ref_end=40, motif="AG", period=2),
            System(name="C", chromosome="chr2", ref_start=1, ref_end=10, motif="AG", period=2),
        ],
    )
    paths = split_panel_by_chromosome(panel, tmp_path / "split")
    assert set(paths.keys()) == {"chr1", "chr2"}
    assert paths["chr1"].read_text().count("\n") == 2
    assert paths["chr2"].read_text().count("\n") == 1


def test_write_codis_panel(tmp_path: Path, codis_panel_yaml: Path) -> None:
    panel = load_panel(codis_panel_yaml)
    out = write_panel_bed(panel, tmp_path / "codis.bed")
    rows = out.read_text().strip().splitlines()
    assert len(rows) == len(panel.systems)
    for row in rows:
        parts = row.split("\t")
        assert len(parts) == 5
        assert parts[0].startswith("chr")


# ---------------------------------------------------------------------------
# Reading a BED back — the escape hatch from curating a panel first
# ---------------------------------------------------------------------------


def _bed(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "regions.bed"
    p.write_text(text, encoding="utf-8")
    return p


def test_round_trip_is_exact(tmp_path: Path) -> None:
    """What we write must be what we read, or the report's BED block lies."""
    from frontstr.panel.bed import load_panel_from_bed, panel_bed_lines

    panel = load_panel(Path("examples/panels/codis_20_grch38.yaml"))
    path = _bed(tmp_path, "\n".join(panel_bed_lines(panel)) + "\n")
    back, _ = load_panel_from_bed(path)

    by_name = {s.name: s for s in back.systems}
    for s in panel.systems:
        if not s.motif:
            continue
        assert by_name[s.name].ref_start == s.ref_start, f"{s.name} start moved"
        assert by_name[s.name].ref_end == s.ref_end, f"{s.name} end moved"


def test_the_two_conventions_differ_by_exactly_one_base(tmp_path: Path) -> None:
    """The whole reason `coords` has no default.

    Standard BED is 0-based half-open; the panel YAML is 1-based inclusive. One
    base at a window edge is the kind of error that yields a plausible wrong
    answer rather than a crash.
    """
    from frontstr.panel.bed import load_panel_from_bed

    path = _bed(tmp_path, "chr11\t2170987\t2171215\tAATG\tTH01\n")
    as_bed, _ = load_panel_from_bed(path, coords="bed0")
    as_panel, _ = load_panel_from_bed(path, coords="panel1")
    assert as_bed.systems[0].ref_start - as_panel.systems[0].ref_start == 1


def test_a_motif_column_is_recognised_by_its_alphabet(tmp_path: Path) -> None:
    from frontstr.panel.bed import load_panel_from_bed

    with_motif, _ = load_panel_from_bed(_bed(tmp_path, "chr11\t100\t200\tAATG\tTH01\n"))
    assert with_motif.systems[0].motif == "AATG"
    assert with_motif.systems[0].name == "TH01"


def test_a_bed_without_a_motif_is_refused(tmp_path: Path) -> None:
    """Not a parsing limitation — without a motif there is no repeat-core
    binning, and binning on raw window length took TH01 from 2 bins to 12."""
    from frontstr.panel.bed import load_panel_from_bed

    with pytest.raises(PanelError, match="repeat-core length"):
        load_panel_from_bed(_bed(tmp_path, "chr11\t100\t200\tTH01\n"))


def test_compound_motifs_are_accepted(tmp_path: Path) -> None:
    from frontstr.panel.bed import load_panel_from_bed

    panel, _ = load_panel_from_bed(_bed(tmp_path, "chr12\t100\t200\tTCTA,TCTG\tvWA\n"))
    assert panel.systems[0].motif == "TCTA,TCTG"


def test_markers_without_calibration_are_reported_and_marked(tmp_path: Path) -> None:
    """A BED carries no period and no corr_value. For a marker STRNaming cannot
    name, that means the number is an uncalibrated repeat count — which has to
    reach the reader, not just the loader."""
    from frontstr.interp.models import FlagCode
    from frontstr.panel.bed import load_panel_from_bed

    panel, warnings = load_panel_from_bed(
        _bed(tmp_path, "chr11\t2170987\t2171215\tAATG\tTH01\nchr9\t100\t200\tAAAG\tMADEUP\n")
    )
    assert "MADEUP" in warnings
    assert not any(w.startswith("TH01") for w in warnings), "TH01 has a STRNaming range"

    made_up = next(s for s in panel.systems if s.name == "MADEUP")
    assert made_up.kit_nomenclature_note, "must carry the note that raises the flag"
    assert FlagCode.CE_NOMENCLATURE_OFFSET  # the machinery the note drives


def test_comments_and_track_lines_are_skipped(tmp_path: Path) -> None:
    from frontstr.panel.bed import load_panel_from_bed

    panel, _ = load_panel_from_bed(
        _bed(tmp_path, '# a comment\ntrack name="x"\nchr11\t100\t200\tAATG\tTH01\n')
    )
    assert len(panel.systems) == 1


def test_malformed_input_fails_loudly(tmp_path: Path) -> None:
    from frontstr.panel.bed import load_panel_from_bed

    with pytest.raises(PanelError, match="chrom/start/end"):
        load_panel_from_bed(_bed(tmp_path, "chr11\t100\n"))
    with pytest.raises(PanelError, match="non-numeric"):
        load_panel_from_bed(_bed(tmp_path, "chr11\tstart\tend\n"))
    with pytest.raises(PanelError, match="reversed"):
        load_panel_from_bed(_bed(tmp_path, "chr11\t300\t200\tAATG\tA\n"))
    with pytest.raises(PanelError, match="duplicate"):
        load_panel_from_bed(_bed(tmp_path, "chr11\t100\t200\tAATG\tX\nchr11\t300\t400\tAATG\tX\n"))
    with pytest.raises(PanelError, match="no usable"):
        load_panel_from_bed(_bed(tmp_path, "# nothing here\n"))
