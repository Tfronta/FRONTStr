"""Tests for the cohort validation harness in ``benchmark/``.

Dev tooling, but the slicing BED it produces decides which reads a two-hour
cohort fetch will ever see, and a region missing from it fails silently rather
than loudly. That is worth a test in the suite CI runs.
"""

from __future__ import annotations

from pathlib import Path

from benchmark.cohort import eligible_basecall, parse_path_list, write_slice_bed

PANEL = Path("examples/panels/codis_20_grch38.yaml")


def _regions(bed: Path) -> dict[str, tuple[str, int, int]]:
    rows = [line.split("\t") for line in bed.read_text().strip().splitlines()]
    return {r[3]: (r[0], int(r[1]), int(r[2])) for r in rows}


def test_slice_bed_includes_the_chry_half_of_amel(tmp_path: Path) -> None:
    """AMEL occupies two regions; dropping the chrY one sex-types males female.

    Regression: the BED was written one row per system, so ``AMEL_Y`` never
    reached samtools. Nothing errored — the slice simply had no chrY reads and
    every male sample called AMEL as ``X``. Caught on HG00112, which reports
    DYS391 and DYS393 (male) but came back ``AMEL=X``.
    """
    bed = write_slice_bed(PANEL, tmp_path / "slice.bed", padding=10_000)
    regions = _regions(bed)

    assert "AMEL" in regions, "the chrX half must be present"
    assert "AMEL_Y" in regions, "the chrY half must be present"
    assert regions["AMEL"][0] == "chrX"
    assert regions["AMEL_Y"][0] == "chrY"


def test_slice_bed_covers_every_panel_system(tmp_path: Path) -> None:
    """Every system gets a region — a missing one is a locus that cannot call."""
    from frontstr.panel.loader import load_panel

    bed = write_slice_bed(PANEL, tmp_path / "slice.bed", padding=10_000)
    regions = _regions(bed)

    for system in load_panel(PANEL).systems:
        assert system.name in regions, f"{system.name} missing from the slicing BED"


def test_slice_bed_pads_both_sides(tmp_path: Path) -> None:
    """The window is widened symmetrically; reads must span repeat plus flanks."""
    padding = 10_000
    bed = write_slice_bed(PANEL, tmp_path / "slice.bed", padding=padding)
    chrom, start, end = _regions(bed)["TPOX"]

    # TPOX: chr2:1489551-1489784 in the panel, 1-based inclusive.
    assert (chrom, start, end) == ("chr2", 1_489_550 - padding, 1_489_784 + padding)


def test_slice_bed_never_starts_before_the_contig(tmp_path: Path) -> None:
    """Padding is clamped at zero rather than producing a negative BED start."""
    bed = write_slice_bed(PANEL, tmp_path / "slice.bed", padding=500_000_000)

    assert all(start == 0 for _, start, _ in _regions(bed).values())


def test_only_r10_dorado_basecalls_are_eligible() -> None:
    """R9 and guppy are rejected: the panel calibration does not describe them."""
    assert eligible_basecall("x/HG00113-ONT-hg38-R10-LSK114-dorado090_sup.phased.bam")
    assert not eligible_basecall("x/GM18501-ONT-hg38-R9-LSK110-guppy-sup-5mC.phased.bam")
    assert not eligible_basecall("x/GM18507-ONT-hg38-R9-LSK110-dorado050_sup.phased.bam")
    assert not eligible_basecall("x/HG00113-ONT-hg38-R10-LSK114-guppy-sup.phased.bam")


def test_path_list_maps_gm_ids_to_the_workbook_spelling(tmp_path: Path) -> None:
    """``GM19038`` in the bucket is ``NA19038`` in the truth workbook."""
    listing = tmp_path / "paths.txt"
    listing.write_text(
        "P/GM19038-ONT-hg38-R10-LSK114-dorado081_sup.phased.bam\n"
        "P/HG00113-ONT-hg38-R10-LSK114-dorado090_sup.phased.bam\n"
        "P/GM18501-ONT-hg38-R9-LSK110-guppy-sup-5mC.phased.bam\n",
        encoding="utf-8",
    )

    by_sample = parse_path_list(listing)

    assert set(by_sample) == {"NA19038", "HG00113"}, "R9/guppy must not appear"
