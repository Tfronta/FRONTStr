"""Tests for the native VCF export.

The file's job is to be read by other people's tools, so the tests read it back
with a real VCF parser rather than asserting on the text we just wrote. A
hand-rolled writer that only its own regex can parse is not an interchange
format.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pysam
import pytest

from frontstr.errors import FrontstrError
from frontstr.exports.vcf import decode_vcf_value, write_run_vcf

CHROM = "chr1"
CONTIG_LEN = 3_000
REF_WINDOW = "ACGT" * 10 + "AGAT" * 8 + "TTCA" * 10  # 40 + 32 + 40 = 112 bp
WINDOW_START = 101  # 1-based
WINDOW_END = WINDOW_START + len(REF_WINDOW) - 1


@pytest.fixture
def reference(tmp_path: Path) -> Path:
    """A tiny indexed FASTA whose window matches :data:`REF_WINDOW`."""
    seq = "N" * (WINDOW_START - 1) + REF_WINDOW
    seq += "N" * (CONTIG_LEN - len(seq))
    fasta = tmp_path / "ref.fa"
    with fasta.open("w") as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")
    pysam.faidx(str(fasta))
    return fasta


def _allele(
    consensus: str,
    n_reads: int,
    *,
    number: float | None = 8.0,
    label: str = "8",
    isfg: str = "[AGAT]8",
    index: int = 0,
    iso: str | None = None,
    hp1: int = 0,
    hp2: int = 0,
    is_deletion: bool = False,
) -> dict[str, Any]:
    return {
        "cluster_index": index,
        "consensus": consensus,
        "length_bp": len(consensus),
        "number": number,
        "number_method": "period_ce",
        "number_label": label,
        "isfg": isfg,
        "is_deletion": is_deletion,
        "n_reads_total": n_reads,
        "n_reads_hp1": hp1,
        "n_reads_hp2": hp2,
        "n_reads_hp_none": n_reads - hp1 - hp2,
        "n_forward": n_reads // 2,
        "n_reverse": n_reads - n_reads // 2,
        "consensus_method": "poa_spoa",
        "iso": {"suffix": iso, "match_type": "exact" if iso else "none"},
        "flags": [],
    }


def _payload(
    called: list[dict[str, Any]],
    *,
    total_reads: int = 40,
    flags: list[dict[str, str]] | None = None,
    marker_type: str = "str",
    name: str = "TEST",
) -> dict[str, Any]:
    return {
        "meta": {"sample_name": "S1", "panel_name": "P", "panel_version": "1"},
        "audit": {"poa_backend": "poa_spoa"},
        "results": [
            {
                "marker_name": name,
                "system": {
                    "chromosome": CHROM,
                    "ref_start": WINDOW_START,
                    "ref_end": WINDOW_END,
                    "motif": "AGAT",
                    "period": 4,
                    "strand": "+",
                    "marker_type": marker_type,
                },
                "total_reads": total_reads,
                "alleles": called,
                "alleles_called": called,
                "flags": flags or [],
            }
        ],
    }


def _read_back(path: Path) -> dict[str, pysam.VariantRecord]:
    with pysam.VariantFile(str(path)) as vf:
        return {rec.id: rec for rec in vf}


# ---------------------------------------------------------------------------
# It has to be a real VCF
# ---------------------------------------------------------------------------


def test_output_parses_as_vcf(tmp_path: Path, reference: Path) -> None:
    out = write_run_vcf(
        _payload([_allele("AGAT" * 20, 30, index=0)]),
        tmp_path / "o.vcf",
        reference_fasta=reference,
    )
    recs = _read_back(out)
    assert list(recs) == ["TEST"]


def test_ref_is_the_actual_reference_sequence(tmp_path: Path, reference: Path) -> None:
    """The claim that makes the file comparable to any other caller's."""
    out = write_run_vcf(
        _payload([_allele("AGAT" * 20, 30)]), tmp_path / "o.vcf", reference_fasta=reference
    )
    rec = _read_back(out)["TEST"]
    with pysam.FastaFile(str(reference)) as fa:
        assert rec.ref == fa.fetch(CHROM, rec.start, rec.stop).upper()
    assert rec.ref == REF_WINDOW


def test_missing_reference_is_refused_not_faked(tmp_path: Path) -> None:
    """A placeholder REF would produce a file that loads and compares wrongly."""
    with pytest.raises(FrontstrError, match="reference FASTA"):
        write_run_vcf(
            _payload([_allele("AGAT" * 20, 30)]), tmp_path / "o.vcf", reference_fasta=None
        )


def test_contig_mismatch_names_the_likely_cause(tmp_path: Path, reference: Path) -> None:
    payload = _payload([_allele("AGAT" * 20, 30)])
    payload["results"][0]["system"]["chromosome"] = "1"
    with pytest.raises(FrontstrError, match="chr prefix"):
        write_run_vcf(payload, tmp_path / "o.vcf", reference_fasta=reference)


def test_records_are_coordinate_sorted(tmp_path: Path, reference: Path) -> None:
    """Panel order is not coordinate order, and an unsorted VCF cannot be indexed."""
    payload = _payload([_allele("AGAT" * 20, 30)], name="LATE")
    late = payload["results"][0]
    early = {**late, "marker_name": "EARLY", "system": {**late["system"]}}
    early["system"]["ref_start"] = 1
    early["system"]["ref_end"] = 50
    payload["results"] = [late, early]

    out = write_run_vcf(payload, tmp_path / "o.vcf", reference_fasta=reference)
    with pysam.VariantFile(str(out)) as vf:
        positions = [r.pos for r in vf]
    assert positions == sorted(positions)


# ---------------------------------------------------------------------------
# Per-allele field alignment — the bug this format is most likely to hide
# ---------------------------------------------------------------------------


def test_per_allele_fields_align_with_ad_when_a_call_equals_ref(
    tmp_path: Path, reference: Path
) -> None:
    """The failure mode: one called allele is identical to the reference.

    Emitted in call order (by depth), MC[0] would describe the *other* allele
    than AD[0], and a benchmark joining those columns would read the locus
    wrong with nothing looking malformed. Every per-allele field is Number=R in
    REF-then-ALT order precisely so that cannot happen.
    """
    ref_call = _allele(REF_WINDOW, 7, label="7", index=1)
    alt_call = _allele(REF_WINDOW + "AGAT", 10, label="9.3", index=0)
    # Call order is by depth, so the ALT (10 reads) comes first.
    out = write_run_vcf(
        _payload([alt_call, ref_call]), tmp_path / "o.vcf", reference_fasta=reference
    )
    rec = _read_back(out)["TEST"]
    sample = rec.samples["S1"]

    assert sample["AD"] == (7, 10), "REF depth first"
    assert sample["MC"] == ("7", "9.3"), "MC must follow the same order as AD"
    assert sample["AL"] == (len(REF_WINDOW), len(REF_WINDOW) + 4)
    assert set(sample["GT"]) == {0, 1}


def test_uncalled_ref_allele_is_missing_not_omitted(tmp_path: Path, reference: Path) -> None:
    """With no call matching REF, index 0 must still hold a slot."""
    out = write_run_vcf(
        _payload(
            [_allele("AGAT" * 20, 18, label="20"), _allele("AGAT" * 19, 16, label="19", index=1)]
        ),
        tmp_path / "o.vcf",
        reference_fasta=reference,
    )
    sample = _read_back(out)["TEST"].samples["S1"]
    assert sample["AD"] == (0, 18, 16)
    assert sample["MC"][0] in (None, "."), "REF slot present but marked missing"
    assert sample["MC"][1:] == ("20", "19")


# ---------------------------------------------------------------------------
# Sequence resolution and encoding
# ---------------------------------------------------------------------------


def test_isoalleles_are_distinct_alts(tmp_path: Path, reference: Path) -> None:
    """Same allele number, different sequence — the whole point of the format."""
    a = _allele("AGAT" * 14 + "AGAC" + "AGAT" * 5, 20, label="20", iso="a")
    b = _allele("AGAT" * 5 + "AGAC" + "AGAT" * 14, 18, label="20", iso="b", index=1)
    out = write_run_vcf(_payload([a, b]), tmp_path / "o.vcf", reference_fasta=reference)
    rec = _read_back(out)["TEST"]
    sample = rec.samples["S1"]

    assert len(rec.alts) == 2, "identical numbers must not collapse into one ALT"
    assert rec.alts[0] != rec.alts[1]
    assert sample["MC"][1:] == ("20", "20")
    assert sample["ISO"][1:] == ("a", "b")


def test_isfg_survives_the_round_trip(tmp_path: Path, reference: Path) -> None:
    """ISFG notation is space-separated; raw spaces are illegal in a VCF field."""
    isfg = "[AATG]6 ATG [AATG]3"
    out = write_run_vcf(
        _payload([_allele("AGAT" * 20, 30, isfg=isfg)]),
        tmp_path / "o.vcf",
        reference_fasta=reference,
    )
    sample = _read_back(out)["TEST"].samples["S1"]
    # Index 0 is REF, which was not called here; the allele is ALT1.
    assert " " not in sample["ISFG"][1]
    assert decode_vcf_value(sample["ISFG"][1]) == isfg


def test_deletion_uses_a_symbolic_alt(tmp_path: Path, reference: Path) -> None:
    out = write_run_vcf(
        _payload([_allele("", 20, is_deletion=True, label="0")]),
        tmp_path / "o.vcf",
        reference_fasta=reference,
    )
    assert _read_back(out)["TEST"].alts == ("<DEL>",)


# ---------------------------------------------------------------------------
# FILTER semantics
# ---------------------------------------------------------------------------


def test_clean_locus_passes(tmp_path: Path, reference: Path) -> None:
    out = write_run_vcf(
        _payload([_allele("AGAT" * 20, 30)]), tmp_path / "o.vcf", reference_fasta=reference
    )
    assert list(_read_back(out)["TEST"].filter) == ["PASS"]


def test_warnings_filter_the_call(tmp_path: Path, reference: Path) -> None:
    flags = [{"code": "low_coverage", "severity": "warn", "message": "x"}]
    out = write_run_vcf(
        _payload([_allele("AGAT" * 20, 12)], flags=flags),
        tmp_path / "o.vcf",
        reference_fasta=reference,
    )
    assert list(_read_back(out)["TEST"].filter) == ["low_coverage"]


def test_informational_flags_do_not_filter(tmp_path: Path, reference: Path) -> None:
    """Resolving an iso-allele is the tool working, not a reason to filter."""
    flags = [{"code": "isoallele", "severity": "info", "message": "x"}]
    out = write_run_vcf(
        _payload([_allele("AGAT" * 20, 30)], flags=flags),
        tmp_path / "o.vcf",
        reference_fasta=reference,
    )
    rec = _read_back(out)["TEST"]
    assert list(rec.filter) == ["PASS"]
    assert rec.info["NOTE"] == ("isoallele",)


# ---------------------------------------------------------------------------
# Exclusions
# ---------------------------------------------------------------------------


def test_amelogenin_is_not_emitted(tmp_path: Path, reference: Path) -> None:
    """X and Y are not nucleotide sequences; a REF/ALT record would be nonsense."""
    payload = _payload([_allele("X", 12, label="X")], marker_type="amel", name="AMEL")
    out = write_run_vcf(payload, tmp_path / "o.vcf", reference_fasta=reference)
    assert _read_back(out) == {}


def test_locus_with_no_reads_is_not_emitted(tmp_path: Path, reference: Path) -> None:
    """An empty record would assert a homozygous-reference call never observed."""
    out = write_run_vcf(_payload([], total_reads=0), tmp_path / "o.vcf", reference_fasta=reference)
    assert _read_back(out) == {}
