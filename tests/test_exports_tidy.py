"""Tests for the cohort-scale tidy dataset.

This is the substrate a concordance study runs on, so the tests care about the
properties that would silently corrupt one: a dropout that vanishes instead of
being recorded, a numeric column typed as string, or run configuration that
cannot be told apart when a cohort spans two calibrations.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest

from frontstr.errors import FrontstrError
from frontstr.exports.tidy import (
    TIDY_COLUMNS,
    build_tidy_rows,
    load_payloads,
    parquet_available,
    write_tidy,
    write_tidy_csv,
    write_tidy_parquet,
)

needs_parquet = pytest.mark.skipif(not parquet_available(), reason="pyarrow not installed")


def _allele(
    label: str,
    number: float | None,
    reads: int,
    *,
    iso: str | None = None,
    absolute: bool = True,
) -> dict[str, Any]:
    return {
        "number": number,
        "number_label": label,
        "number_method": "period_ce" if absolute else "none",
        "number_is_absolute": absolute,
        "isfg": f"[AGAT]{label}",
        "iso": {"suffix": iso, "match_type": "exact" if iso else "none"},
        "length_bp": 200,
        "n_reads_total": reads,
        "n_reads_hp1": 0,
        "n_reads_hp2": 0,
        "n_reads_absorbed": 0,
        "fraction": 0.5,
        "consensus_method": "poa_spoa",
        "consensus": "AGAT" * 50,
    }


def _marker(
    name: str,
    called: list[dict[str, Any]],
    *,
    call_rule: str = "heterozygous",
    flags: list[dict[str, str]] | None = None,
    total_reads: int = 40,
) -> dict[str, Any]:
    return {
        "marker_name": name,
        "call_rule": call_rule,
        "total_reads": total_reads,
        "alleles_called": called,
        "flags": flags or [],
    }


def _payload(
    sample: str,
    markers: list[dict[str, Any]],
    *,
    panel_version: str = "1.0",
    stutter_model: str = "m1",
) -> dict[str, Any]:
    return {
        "meta": {
            "sample_name": sample,
            "panel_name": "CODIS",
            "panel_version": panel_version,
            "run_id": "R1",
        },
        "audit": {"poa_backend": "poa_spoa", "stutter_model_version": stutter_model},
        "results": markers,
    }


# ---------------------------------------------------------------------------
# Shape
# ---------------------------------------------------------------------------


def test_one_row_per_sample_marker_allele() -> None:
    payload = _payload(
        "S1",
        [
            _marker("M1", [_allele("8", 8.0, 20), _allele("9", 9.0, 18)]),
            _marker("M2", [_allele("11", 11.0, 30)], call_rule="homozygous"),
        ],
    )
    rows = build_tidy_rows([payload])
    assert len(rows) == 3
    assert [(r["marker"], r["allele_index"]) for r in rows] == [
        ("M1", 1),
        ("M1", 2),
        ("M2", 1),
    ]


def test_dropout_is_a_row_not_an_absence() -> None:
    """Omitting it makes a dropout indistinguishable from a marker never run."""
    payload = _payload("S1", [_marker("M1", [], call_rule="no_data", total_reads=0)])
    rows = build_tidy_rows([payload])
    assert len(rows) == 1
    row = rows[0]
    assert row["called"] is False
    assert row["allele_index"] == 0
    assert row["allele_number"] is None
    assert row["call_rule"] == "no_data"
    assert row["marker"] == "M1"


def test_multiple_payloads_combine_into_one_table() -> None:
    rows = build_tidy_rows(
        [
            _payload("S1", [_marker("M1", [_allele("8", 8.0, 20)])]),
            _payload("S2", [_marker("M1", [_allele("9", 9.0, 22)])]),
        ]
    )
    assert [r["sample"] for r in rows] == ["S1", "S2"]


def test_every_declared_column_is_populated() -> None:
    rows = build_tidy_rows([_payload("S1", [_marker("M1", [_allele("8", 8.0, 20)])])])
    assert set(rows[0]) == set(TIDY_COLUMNS)


# ---------------------------------------------------------------------------
# Content a concordance study depends on
# ---------------------------------------------------------------------------


def test_allele_number_is_carried_as_a_number() -> None:
    rows = build_tidy_rows([_payload("S1", [_marker("M1", [_allele("9.3", 9.3, 20)])])])
    assert rows[0]["allele_number"] == 9.3
    assert rows[0]["allele_label"] == "9.3"


def test_non_comparable_numbers_are_marked() -> None:
    """AMEL's X/Y have no allele number; joining them to a CE table is an error."""
    rows = build_tidy_rows(
        [_payload("S1", [_marker("AMEL", [_allele("X", None, 12, absolute=False)])])]
    )
    assert rows[0]["number_is_absolute"] is False
    assert rows[0]["allele_number"] is None
    assert rows[0]["allele_label"] == "X"


def test_marker_flags_and_review_status_travel_with_every_allele() -> None:
    payload = _payload(
        "S1",
        [
            _marker(
                "M1",
                [_allele("8", 8.0, 20), _allele("9", 9.0, 18)],
                flags=[{"code": "low_coverage", "severity": "warn", "message": "x"}],
            )
        ],
    )
    rows = build_tidy_rows([payload])
    assert all(r["marker_flags"] == "low_coverage" for r in rows)
    assert all(r["needs_review"] for r in rows)


def test_informational_flags_do_not_mark_a_row_for_review() -> None:
    payload = _payload(
        "S1",
        [
            _marker(
                "M1",
                [_allele("8", 8.0, 20)],
                flags=[{"code": "isoallele", "severity": "info", "message": "x"}],
            )
        ],
    )
    rows = build_tidy_rows([payload])
    assert rows[0]["needs_review"] is False
    assert rows[0]["marker_flags"] == "isoallele"


def test_run_configuration_distinguishes_cohort_halves() -> None:
    """A 150-sample cohort is not collected in an afternoon.

    When half was called under one stutter model, the dataset has to say so.
    """
    rows = build_tidy_rows(
        [
            _payload("S1", [_marker("M1", [_allele("8", 8.0, 20)])], stutter_model="old"),
            _payload("S2", [_marker("M1", [_allele("8", 8.0, 20)])], stutter_model="new"),
        ]
    )
    assert [r["stutter_model"] for r in rows] == ["old", "new"]
    assert all(r["poa_backend"] == "poa_spoa" for r in rows)


def test_isoallele_suffix_is_carried() -> None:
    rows = build_tidy_rows([_payload("S1", [_marker("M1", [_allele("14", 14.0, 20, iso="b")])])])
    assert rows[0]["iso_suffix"] == "b"
    assert rows[0]["iso_match"] == "exact"


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------


def test_csv_has_the_declared_header_in_order(tmp_path: Path) -> None:
    rows = build_tidy_rows([_payload("S1", [_marker("M1", [_allele("8", 8.0, 20)])])])
    out = write_tidy_csv(rows, tmp_path / "t.csv")
    with out.open() as fh:
        assert next(csv.reader(fh)) == list(TIDY_COLUMNS)


def test_consensus_is_the_last_column() -> None:
    """It is 250 bp wide; anywhere else it pushes the useful columns off screen."""
    assert TIDY_COLUMNS[-1] == "consensus"


# ---------------------------------------------------------------------------
# Parquet
# ---------------------------------------------------------------------------


@needs_parquet
def test_parquet_types_are_explicit_not_inferred(tmp_path: Path) -> None:
    """The failure this prevents: allele_number typed as string, so "10" < "9".

    Inference reads the first rows; a marker that dropped out puts nulls there.
    """
    import pyarrow.parquet as pq

    rows = build_tidy_rows(
        [
            _payload("S1", [_marker("M0", [], call_rule="no_data", total_reads=0)]),
            _payload("S2", [_marker("M1", [_allele("9.3", 9.3, 20)])]),
        ]
    )
    out = write_tidy_parquet(rows, tmp_path / "t.parquet")
    schema = pq.read_schema(out)
    assert schema.field("allele_number").type == "double"
    assert schema.field("reads").type == "int32"
    assert schema.field("called").type == "bool"
    assert schema.field("sample").type == "string"


@needs_parquet
def test_parquet_round_trips_values(tmp_path: Path) -> None:
    import pyarrow.parquet as pq

    rows = build_tidy_rows(
        [
            _payload("S1", [_marker("M1", [_allele("9.3", 9.3, 20)])]),
            _payload("S2", [_marker("M1", [], call_rule="no_data", total_reads=0)]),
        ]
    )
    table = pq.read_table(write_tidy_parquet(rows, tmp_path / "t.parquet"))
    data = table.to_pydict()
    assert data["allele_number"] == [9.3, None]
    assert data["called"] == [True, False]
    assert data["sample"] == ["S1", "S2"]


@needs_parquet
def test_write_tidy_emits_both_formats(tmp_path: Path) -> None:
    written = write_tidy([_payload("S1", [_marker("M1", [_allele("8", 8.0, 20)])])], tmp_path)
    assert [p.suffix for p in written] == [".csv", ".parquet"]
    assert all(p.exists() for p in written)


def test_write_tidy_degrades_to_csv_without_pyarrow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A machine without pyarrow should still get a usable dataset."""
    monkeypatch.setattr("frontstr.exports.tidy.parquet_available", lambda: False)
    written = write_tidy([_payload("S1", [_marker("M1", [_allele("8", 8.0, 20)])])], tmp_path)
    assert [p.suffix for p in written] == [".csv"]


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def test_load_payloads_reads_run_json(tmp_path: Path) -> None:
    path = tmp_path / "run.json"
    path.write_text(json.dumps(_payload("S1", [_marker("M1", [_allele("8", 8.0, 20)])])))
    assert [p["meta"]["sample_name"] for p in load_payloads([path])] == ["S1"]


def test_load_payloads_rejects_a_json_that_is_not_a_run(tmp_path: Path) -> None:
    """Pointing --from-dir at a directory of unrelated JSON must say so."""
    path = tmp_path / "other.json"
    path.write_text('{"something": "else"}')
    with pytest.raises(FrontstrError, match="not a FRONTStr run JSON"):
        list(load_payloads([path]))


def test_load_payloads_reports_unreadable_files(tmp_path: Path) -> None:
    path = tmp_path / "broken.json"
    path.write_text("{not json")
    with pytest.raises(FrontstrError, match="Cannot read run JSON"):
        list(load_payloads([path]))
