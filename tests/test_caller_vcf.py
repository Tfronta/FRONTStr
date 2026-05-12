"""Tests for cyvcf2-based LongTR VCF parsing.

These tests construct hand-crafted minimal VCF files to exercise every branch
of the parser without requiring an actual LongTR install.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from frontstr.caller.vcf import (
    LongTRResult,
    parse_longtr_vcf,
)
from frontstr.caller.vcf import _parse_allreads as parse_allreads
from frontstr.caller.vcf import _parse_pdp as parse_pdp
from frontstr.errors import CallerError

_VCF_HEADER = """##fileformat=VCFv4.2
##INFO=<ID=MOTIF,Number=.,Type=String,Description="Motif">
##INFO=<ID=PERIOD,Number=.,Type=String,Description="Period">
##INFO=<ID=BPDIFFS,Number=A,Type=Integer,Description="bp diff per ALT">
##INFO=<ID=INEXACT_ALLELE,Number=A,Type=Integer,Description="POA-derived ALT flag">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=Q,Number=1,Type=Float,Description="Posterior probability">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=PDP,Number=1,Type=String,Description="HP1|HP2 counts">
##FORMAT=<ID=ALLREADS,Number=1,Type=String,Description="bpdiff|reads pairs">
##contig=<ID=chr1>
##contig=<ID=chr2>
##contig=<ID=chr11>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def _write_vcf(path: Path, records: list[str]) -> Path:
    path.write_text(_VCF_HEADER + "\n".join(records) + "\n")
    return path


def test_parse_simple_heterozygous(tmp_path: Path) -> None:
    rec = (
        "chr11\t2171100\tTH01\tAATG\tAATGAATG,AATGAATGAATG\t."
        "\tPASS"
        "\tMOTIF=AATG;PERIOD=4;BPDIFFS=4,8;INEXACT_ALLELE=0,0"
        "\tGT:Q:DP:PDP:ALLREADS"
        "\t1/2:0.99:120:60|55:0|3;4|62;8|55"
    )
    vcf = _write_vcf(tmp_path / "th01.vcf", [rec])
    results = parse_longtr_vcf(vcf)
    assert len(results) == 1
    r: LongTRResult = results[0]
    assert r.marker_name == "TH01"
    assert r.chrom == "chr11"
    assert r.motif == "AATG"
    assert r.period == "4"
    assert len(r.alleles) == 3
    assert r.alleles[0].sequence == "AATG"
    assert r.alleles[0].bp_diff == 0
    assert r.alleles[1].bp_diff == 4
    assert r.alleles[2].bp_diff == 8
    assert not any(a.inexact for a in r.alleles)
    assert not any(a.is_deletion for a in r.alleles)

    call = r.samples["S1"]
    assert call.gt_indices == (1, 2)
    assert call.phased is False
    assert call.posterior == pytest.approx(0.99)
    assert call.depth == 120
    assert call.pdp_hp1 == 60
    assert call.pdp_hp2 == 55
    assert call.allreads == {0: 3, 4: 62, 8: 55}


def test_parse_deletion_allele(tmp_path: Path) -> None:
    rec = (
        "chr1\t1000\tFGA\tAGAT\t<DEL>\t."
        "\tPASS"
        "\tMOTIF=AGAT;PERIOD=4;BPDIFFS=-4;INEXACT_ALLELE=0"
        "\tGT:Q:DP:PDP:ALLREADS"
        "\t0/1:0.91:30:15|15:0|18;-4|12"
    )
    vcf = _write_vcf(tmp_path / "fga.vcf", [rec])
    r = parse_longtr_vcf(vcf)[0]
    assert len(r.alleles) == 2
    assert r.alleles[1].sequence == "<DEL>"
    assert r.alleles[1].is_deletion is True
    assert r.alleles[1].bp_diff == -4


def test_parse_inexact_allele_flag(tmp_path: Path) -> None:
    rec = (
        "chr1\t500\tVNTR1\tACGTACGTACGT\tACGTACGTACGTACGT\t."
        "\tPASS"
        "\tMOTIF=ACGT;PERIOD=4;BPDIFFS=4;INEXACT_ALLELE=1"
        "\tGT:Q:DP:PDP:ALLREADS"
        "\t0/1:0.95:50:25|24:0|24;4|20"
    )
    vcf = _write_vcf(tmp_path / "vntr.vcf", [rec])
    r = parse_longtr_vcf(vcf)[0]
    assert r.alleles[0].inexact is False
    assert r.alleles[1].inexact is True


def test_parse_phased_genotype(tmp_path: Path) -> None:
    rec = (
        "chr1\t1000\tM\tAGAT\tAGATAGAT\t."
        "\tPASS"
        "\tMOTIF=AGAT;PERIOD=4;BPDIFFS=4;INEXACT_ALLELE=0"
        "\tGT:Q:DP:PDP:ALLREADS"
        "\t0|1:0.98:100:50|45:0|55;4|45"
    )
    vcf = _write_vcf(tmp_path / "phased.vcf", [rec])
    r = parse_longtr_vcf(vcf)[0]
    call = r.samples["S1"]
    assert call.phased is True
    assert call.gt_indices == (0, 1)


def test_parse_missing_genotype(tmp_path: Path) -> None:
    rec = (
        "chr1\t1000\tM\tAGAT\tAGATAGAT\t."
        "\tPASS"
        "\tMOTIF=AGAT;PERIOD=4;BPDIFFS=4;INEXACT_ALLELE=0"
        "\tGT:Q:DP:PDP:ALLREADS"
        "\t./.:.:.:.:."
    )
    vcf = _write_vcf(tmp_path / "miss.vcf", [rec])
    r = parse_longtr_vcf(vcf)[0]
    call = r.samples["S1"]
    assert call.gt_indices is None
    assert call.posterior is None
    assert call.depth == 0
    assert call.pdp_hp1 == 0
    assert call.pdp_hp2 == 0
    assert call.allreads == {}


def test_parse_multi_motif_locus(tmp_path: Path) -> None:
    rec = (
        "chr3\t45582240\tD3S1358\tTCTA\tTCTATCTA\t."
        "\tPASS"
        "\tMOTIF=TCTA,TCTG;PERIOD=4,4;BPDIFFS=4;INEXACT_ALLELE=0"
        "\tGT:Q:DP:PDP:ALLREADS"
        "\t0/1:0.99:100:50|50:0|50;4|50"
    )
    vcf = _write_vcf(tmp_path / "d3.vcf", [rec])
    r = parse_longtr_vcf(vcf)[0]
    assert r.motif == "TCTA,TCTG"
    assert r.period == "4,4"


def test_parse_no_alt(tmp_path: Path) -> None:
    """Locus where every read supports the reference allele (no ALT alleles)."""
    rec = (
        "chr1\t1000\tHOM\tAGAT\t.\t."
        "\tPASS"
        "\tMOTIF=AGAT;PERIOD=4"
        "\tGT:Q:DP:PDP:ALLREADS"
        "\t0/0:0.99:50:25|25:0|50"
    )
    vcf = _write_vcf(tmp_path / "hom.vcf", [rec])
    r = parse_longtr_vcf(vcf)[0]
    assert len(r.alleles) == 1
    assert r.alleles[0].sequence == "AGAT"
    assert r.alleles[0].bp_diff == 0


def test_parse_missing_file() -> None:
    with pytest.raises(CallerError, match="VCF not found"):
        parse_longtr_vcf(Path("/nonexistent/file.vcf"))


def test_parse_allreads_helper() -> None:
    assert parse_allreads("0|45;-3|2;6|18") == {0: 45, -3: 2, 6: 18}
    assert parse_allreads("") == {}
    assert parse_allreads(".") == {0: 0} or parse_allreads(".") == {}  # accept either
    assert parse_allreads("garbage") == {}


def test_parse_pdp_helper() -> None:
    assert parse_pdp("12|28") == (12, 28)
    assert parse_pdp("") == (0, 0)
    assert parse_pdp(".") == (0, 0)
    assert parse_pdp("garbage") == (0, 0)
