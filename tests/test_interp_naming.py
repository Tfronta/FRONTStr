"""Unit tests for STRNaming-backed allele naming.

These run off the committed slice cache and need neither the ONT BAMs nor a
reference FASTA — that hermeticity is itself part of what is under test.
"""

from __future__ import annotations

import pytest

from frontstr.interp.naming import (
    CACHE_PATH,
    NameStatus,
    StrNamer,
    _extract_range,
    _strip_dropped_options,
    default_namer,
    reverse_complement,
)


@pytest.fixture(scope="module")
def namer() -> StrNamer:
    return StrNamer.from_cache()


def _cache_rows() -> list[list[str]]:
    return [
        line.split("\t")
        for line in CACHE_PATH.read_text(encoding="utf-8").splitlines()[1:]
        if line.strip()
    ]


def _reference_window(row: list[str], pad: int = 120) -> str:
    """The marker's reference sequence with flank, as a cluster consensus would be."""
    _, _, start, end, _, _, slice_start, slice_seq = row
    lo, hi = min(int(start), int(end)), max(int(start), int(end))
    offset = int(slice_start)
    return slice_seq[max(0, lo - offset - pad) : hi - offset + 1 + pad]


def test_reverse_complement_round_trips() -> None:
    assert reverse_complement("AACGT") == "ACGTT"
    assert reverse_complement(reverse_complement("ACGTACGGT")) == "ACGTACGGT"


def test_cache_is_committed_and_non_empty() -> None:
    rows = _cache_rows()
    assert len(rows) >= 20, "the bundled slice cache looks truncated"
    assert all(len(r) == 8 for r in rows)


@pytest.mark.parametrize("row", _cache_rows(), ids=lambda r: r[0])
def test_reference_allele_names_match_ref_ce(namer: StrNamer, row: list[str]) -> None:
    """Naming the reference window must reproduce STRNaming's own ``ref_ce``.

    This is the calibration anchor for the whole module: the reference allele's
    designation is published, so any drift in extraction, orientation or store
    construction shows up here rather than in a sample profile.
    """
    marker, ref_ce = row[0], row[4]
    result = namer.name(marker, _reference_window(row))
    assert result.status is NameStatus.OK, f"{marker}: {result.status}"
    assert result.ce == float(ref_ce)
    assert result.name.startswith(f"CE{ref_ce}_")


def test_minus_strand_marker_is_oriented(namer: StrNamer) -> None:
    """vWA's range is minus-strand; naming must not silently use the forward read."""
    row = next(r for r in _cache_rows() if r[0] == "vWA")
    result = namer.name("vWA", _reference_window(row))
    # Forward-strand motifs per ISFG Recommendation 1, not the panel's TCTA/TCTG.
    assert "AGAT" in result.name
    assert result.ce == 17.0


def test_unknown_marker_reports_no_range(namer: StrNamer) -> None:
    assert namer.has_range("DYS393") is False
    result = namer.name("DYS393", "ACGT" * 40)
    assert result.status is NameStatus.NO_RANGE
    assert result.ce is None


def test_empty_consensus_reports_empty(namer: StrNamer) -> None:
    assert namer.name("vWA", "").status is NameStatus.EMPTY


def test_consensus_without_the_range_is_rejected(namer: StrNamer) -> None:
    """A sloppy anchor hit must not become a confident-looking wrong name."""
    result = namer.name("vWA", "N" * 300)
    assert result.status is NameStatus.ANCHOR_NOT_FOUND
    assert result.ce is None


def test_naming_never_raises_on_garbage(namer: StrNamer) -> None:
    for junk in ("", "A", "ACGT", "N" * 5000, "acgt" * 50):
        assert namer.name("D21S11", junk).ce is None or True  # must simply not raise


def test_default_namer_is_cached() -> None:
    assert default_namer() is default_namer()


class TestExtractRange:
    """Anchored extraction — the part that coordinate slicing got wrong."""

    def test_extracts_between_anchors(self) -> None:
        left, right = "ACGTACGTAC", "TTGGTTGGTT"
        assert _extract_range(f"CCC{left}PAYLOAD{right}GGG", left, right) == "PAYLOAD"

    def test_tolerates_sequencing_error_in_the_anchor(self) -> None:
        left, right = "ACGTACGTACGTACGTACGTACGTACGTAC", "TTGGTTGGTTGGTTGGTTGGTTGGTTGGTT"
        noisy_left = left[:10] + "A" + left[11:]  # one substitution
        assert _extract_range(f"{noisy_left}PAYLOAD{right}", left, right) == "PAYLOAD"

    def test_rejects_when_an_anchor_is_absent(self) -> None:
        assert _extract_range("N" * 200, "ACGTACGTACGTACGTACGTACGTACGTAC", "TTGG") is None

    def test_rejects_crossed_anchors(self) -> None:
        left, right = "ACGTACGTAC", "TTGGTTGGTT"
        assert _extract_range(f"{right}{left}", left, right) is None

    def test_captures_an_insertion_the_reference_does_not_have(self) -> None:
        """The TPOX failure mode: extra repeat units must land inside the range."""
        left, right = "ACGTACGTACGTACGTACGTACGTACGTAC", "TTGGTTGGTTGGTTGGTTGGTTGGTTGGTT"
        expanded = "AATG" * 11
        assert _extract_range(f"{left}{expanded}{right}", left, right) == expanded


def test_strip_dropped_options_removes_only_the_uas_limit() -> None:
    assert _strip_dropped_options("limit=168") == ""
    assert _strip_dropped_options("") == ""
    assert _strip_dropped_options("limit=197,foo=1") == "foo=1"
    assert _strip_dropped_options("foo=1") == "foo=1"


def test_uas_limit_guard_is_disabled(namer: StrNamer) -> None:
    """DXS8378's range is exactly its UAS truncation limit.

    Left in place, the guard makes the marker unnameable even on the reference
    (``CE?_TODO_UAS_INCOMPLETE_SEQUENCE``). FRONTStr never truncates, so the
    guard is a pure false positive here.
    """
    row = next(r for r in _cache_rows() if r[0] == "DXS8378")
    result = namer.name("DXS8378", _reference_window(row))
    assert result.status is NameStatus.OK
    assert result.ce == 10.0
