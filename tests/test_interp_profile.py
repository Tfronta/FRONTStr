"""End-to-end tests for the interp orchestrator using synthetic BAMs."""

from __future__ import annotations

from pathlib import Path

from frontstr.evidence.cluster import cluster_observations
from frontstr.evidence.pileup import pileup_locus
from frontstr.interp.models import AlleleStatus, CallRule, TriType
from frontstr.interp.profile import interpret_marker, interpret_run
from frontstr.panel.models import Panel, System

from .conftest import SYNTH_CHROM, SYNTH_TR_END, SYNTH_TR_START, SynthRead, write_synth_bam


def _synth_system() -> System:
    return System(
        name="SYNTH",
        chromosome=SYNTH_CHROM,
        ref_start=SYNTH_TR_START + 1,
        ref_end=SYNTH_TR_END,
        motif="AGAT",
        period=4,
    )


def test_interpret_marker_heterozygous_high_coverage(tmp_path: Path) -> None:
    """30x CE12 + 25x CE11 + 3x CE10 (real-world forensic-grade coverage) → het + CE10 stutter."""
    specs = (
        [SynthRead(name=f"a{i}", n_repeats=12, hp=1) for i in range(30)]
        + [SynthRead(name=f"b{i}", n_repeats=11, hp=2) for i in range(25)]
        + [SynthRead(name=f"c{i}", n_repeats=10) for i in range(3)]
    )
    bam = write_synth_bam(tmp_path / "het30x.bam", specs)
    obs = pileup_locus(bam, SYNTH_CHROM, SYNTH_TR_START, SYNTH_TR_END, min_mapq=20)
    clusters = cluster_observations(obs)
    result = interpret_marker(system=_synth_system(), clusters=clusters)
    assert result.total_reads == 58
    assert result.call_rule == CallRule.HETEROZYGOUS
    assert result.tri_type == TriType.NONE
    called_ce = sorted(a.ce for a in result.alleles_called if a.ce is not None)
    assert called_ce == [11.0, 12.0]
    statuses = {round(a.ce or 0, 1): a.status for a in result.alleles}
    # CE10 is below the calling threshold (~5%) — either stutter or artefact is forensically valid
    assert statuses[10.0] in {AlleleStatus.STUTTER, AlleleStatus.ARTEFACT}


def test_interpret_marker_low_coverage_phantom_not_mixture(
    synth_bam_heterozygous: Path,
) -> None:
    """Bug #6 fix: at 11x total, the 2-read CE10 cluster (18%) clears the
    fractional calling threshold but falls below the absolute ``min_reads_third``
    floor (default 5). It is an ONT basecaller phantom, not a 3rd contributor,
    so FRONTStr calls heterozygous instead of raising a false mixture flag."""
    obs = pileup_locus(
        synth_bam_heterozygous,
        SYNTH_CHROM,
        SYNTH_TR_START,
        SYNTH_TR_END,
        min_mapq=20,
    )
    clusters = cluster_observations(obs)
    result = interpret_marker(system=_synth_system(), clusters=clusters)
    assert result.total_reads == 11
    assert result.call_rule == CallRule.HETEROZYGOUS
    assert result.tri_type == TriType.NONE


def test_interpret_marker_homozygous(synth_bam_homozygous: Path) -> None:
    obs = pileup_locus(
        synth_bam_homozygous,
        SYNTH_CHROM,
        SYNTH_TR_START,
        SYNTH_TR_END,
        min_mapq=20,
    )
    clusters = cluster_observations(obs)
    result = interpret_marker(system=_synth_system(), clusters=clusters)
    assert result.call_rule == CallRule.HOMOZYGOUS
    assert len(result.alleles_called) == 1
    assert result.alleles_called[0].ce == 12.0
    assert result.alleles_called[0].allele_numeric == 12.0
    assert result.alleles_called[0].allele_numeric_source == "period_ce"
    # Strand balance preserved in the called allele
    assert result.alleles_called[0].n_forward + result.alleles_called[0].n_reverse == 8


def test_interpret_marker_low_mapq_filtered(synth_bam_low_mapq: Path) -> None:
    """Low-MAPQ reads must be filtered before clustering — only CE12 reaches call."""
    obs = pileup_locus(
        synth_bam_low_mapq,
        SYNTH_CHROM,
        SYNTH_TR_START,
        SYNTH_TR_END,
        min_mapq=20,
    )
    clusters = cluster_observations(obs)
    result = interpret_marker(system=_synth_system(), clusters=clusters)
    assert result.total_reads == 3
    assert result.call_rule == CallRule.HOMOZYGOUS
    assert result.alleles_called[0].ce == 12.0


def test_interpret_run_empty_locus_returns_no_data(tmp_path: Path) -> None:
    """A marker with zero reads must produce a NO_DATA result, not raise."""
    bam = write_synth_bam(tmp_path / "empty.bam", [])
    panel = Panel(
        name="t",
        version="0",
        systems=[_synth_system()],
    )
    results = interpret_run(bam=bam, panel=panel)
    assert len(results) == 1
    assert results[0].call_rule == CallRule.NO_DATA
    assert results[0].total_reads == 0


def test_interpret_marker_isfg_compression(synth_bam_homozygous: Path) -> None:
    """ISFG nomenclature should compress 12 AGATs to ``[AGAT]12``."""
    obs = pileup_locus(
        synth_bam_homozygous,
        SYNTH_CHROM,
        SYNTH_TR_START,
        SYNTH_TR_END,
        min_mapq=20,
    )
    clusters = cluster_observations(obs)
    result = interpret_marker(system=_synth_system(), clusters=clusters)
    assert result.alleles_called[0].isfg == "[AGAT]12"


def test_interpret_marker_uses_longtr_inexact_flag(synth_bam_heterozygous: Path) -> None:
    """When LongTR flags the called allele as INEXACT_ALLELE, status is upgraded."""
    from frontstr.caller.vcf import LongTRAlleleSpec, LongTRResult, LongTRSampleCall

    obs = pileup_locus(
        synth_bam_heterozygous,
        SYNTH_CHROM,
        SYNTH_TR_START,
        SYNTH_TR_END,
        min_mapq=20,
    )
    clusters = cluster_observations(obs)
    ref = "AGAT" * 12
    longtr = LongTRResult(
        marker_name="SYNTH",
        chrom=SYNTH_CHROM,
        pos=SYNTH_TR_START + 1,
        motif="AGAT",
        period="4",
        alleles=[
            LongTRAlleleSpec(sequence=ref, bp_diff=0, inexact=False, is_deletion=False),
            LongTRAlleleSpec(sequence="AGAT" * 11, bp_diff=-4, inexact=True, is_deletion=False),
        ],
        samples={
            "S1": LongTRSampleCall(
                sample="S1",
                gt_indices=(0, 1),
                phased=False,
                posterior=0.99,
                depth=9,
            )
        },
    )
    result = interpret_marker(
        system=_synth_system(),
        clusters=clusters,
        longtr=longtr,
    )
    inexact_alleles = [a for a in result.alleles if a.longtr_inexact]
    assert any(a.consensus == "AGAT" * 11 for a in inexact_alleles)
