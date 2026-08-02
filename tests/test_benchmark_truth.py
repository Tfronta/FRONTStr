"""Tests for reading the external truth workbook.

The workbook is a merge of several pipelines and does not use one notation
throughout. HipSTR and longTR report an allele in bases, which the sheet
divides by the motif length, so three extra bases on a tetranucleotide are
stored as ``.75`` where ISFG writes ``.3``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from benchmark.truth import (
    MOTIF_LENGTHS,
    canonical_allele,
    canonical_genotype,
    isfg_from_repeat_units,
)
from frontstr.panel.loader import load_panel

PANEL = Path("examples/panels/codis_20_grch38.yaml")


@pytest.mark.parametrize(
    ("stored", "motif_length", "expected"),
    [
        (18.75, 4, "18.3"),  # 75 bp core: 18 repeats and 3 extra bases
        (9.25, 4, "9.1"),
        (10.5, 4, "10.2"),
        (11.6667, 3, "11.2"),
        (11.3333, 3, "11.1"),
    ],
)
def test_repeat_units_become_isfg(stored: float, motif_length: int, expected: str) -> None:
    assert isfg_from_repeat_units(stored, motif_length) == expected


@pytest.mark.parametrize("stored", [9.1, 9.2, 9.3, 12.0, 14.0])
def test_an_isfg_cell_is_left_alone(stored: float) -> None:
    """.1/.2/.3 cannot come out of a division by four, so they are ISFG."""
    assert isfg_from_repeat_units(stored, 4) is None


def test_a_pentanucleotide_is_never_reinterpreted() -> None:
    """At five, ``.2`` is valid in both notations and means different alleles.

    Two extra bases in ISFG, one in repeat units. Nothing distinguishes them,
    so the cell has to stand as written.
    """
    assert isfg_from_repeat_units(9.2, 5) is None
    assert isfg_from_repeat_units(9.4, 5) is None


def test_an_unrecognised_suffix_is_left_alone() -> None:
    """.45 is neither notation. Reporting it beats inventing an allele."""
    assert isfg_from_repeat_units(9.45, 4) is None


def test_the_workbook_cell_that_started_this() -> None:
    """HG00119 vWA: the sheet stores 18.75 and Excel displays it as 18.8.

    Rounding it to one decimal scored FRONTStr's correct 18.3 as a mismatch.
    """
    assert canonical_allele(18.75, motif_length=4) == "18.3"
    assert canonical_allele(18.75) == "18.8", "unconverted, this is the old reading"


def test_a_stored_quarter_is_not_rounded_into_a_real_allele() -> None:
    """9.25 renders as "9.2" under round-half-even, and 9.2 is a callable
    allele. That is the dangerous direction: it can manufacture agreement."""
    assert canonical_allele(9.25) == "9.2", "the rounding this replaces"
    assert canonical_allele(9.25, motif_length=4) == "9.1"


def test_frontstr_output_is_never_reinterpreted() -> None:
    """FRONTStr emits ISFG, so it is read with no motif length at all."""
    assert canonical_genotype([9.1, 12]) == ("9.1", "12")


def test_conversions_are_reported() -> None:
    converted: list[tuple[float, str]] = []
    canonical_genotype([18.75, 14], motif_length=4, converted=converted)
    assert converted == [(18.75, "18.3")]


def test_motif_lengths_match_the_shipped_panel() -> None:
    """The map is a copy of the panel's motifs, so it can drift. This is what
    catches it."""
    panel = load_panel(PANEL)
    systems = getattr(panel, "markers", None) or panel.systems
    for system in systems:
        if system.name not in MOTIF_LENGTHS:
            continue
        lengths = {len(m) for m in system.motif.split(",") if m}
        assert lengths == {MOTIF_LENGTHS[system.name]}, system.name
