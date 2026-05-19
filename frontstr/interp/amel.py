"""Sex typing from the AMEL (Amelogenin) locus.

AMEL is not an STR. Forensic sex determination works by counting reads that
align to the AMELX region (chrX) and the AMELY region (chrY) separately:

- Reads only at chrX  → female (X)
- Reads at chrX + chrY → male (X/Y)
- No reads            → no_data

The X and Y alleles are reported as ISFG-style "X" / "Y" strings; CE is not
applicable (``None``). The AMELY coordinates for the System must be set via
the ``y_chromosome``, ``y_ref_start``, and ``y_ref_end`` fields.
"""

from __future__ import annotations

from pathlib import Path

from frontstr.evidence.pileup import pileup_locus
from frontstr.interp.models import Allele, AlleleStatus, CallRule, MarkerResult, TriType
from frontstr.panel.models import System

_DEFAULT_MIN_READS = 3


def interpret_amel(
    system: System,
    bam: Path,
    *,
    min_mapq: int = 20,
    min_reads: int = _DEFAULT_MIN_READS,
    reference_fasta: Path | None = None,
) -> MarkerResult:
    """Sex-type one sample from AMELX and AMELY read counts.

    Args:
        system: The AMEL System. Should have ``y_chromosome``, ``y_ref_start``,
            and ``y_ref_end`` set for Y pileup. Without them only chrX is
            checked, producing a female-only call.
        bam: Indexed sample BAM or CRAM aligned to a GRCh38 reference that
            includes both chrX and chrY (pseudo-autosomal regions masked).
        min_mapq: Minimum mapping quality for pileup reads.
        min_reads: Minimum reads on a chromosome to call that allele present.
        reference_fasta: Reference FASTA path; required when ``bam`` is a CRAM.
    """
    x_count = len(
        _safe_pileup(bam, system.chromosome, system.ref_start - 1, system.ref_end, min_mapq,
                     reference_fasta=reference_fasta)
    )
    y_count = 0
    if system.y_chromosome and system.y_ref_start and system.y_ref_end:
        y_count = len(
            _safe_pileup(
                bam, system.y_chromosome, system.y_ref_start - 1, system.y_ref_end, min_mapq,
                reference_fasta=reference_fasta,
            )
        )

    total_reads = x_count + y_count
    if total_reads == 0:
        return _no_data(system)

    alleles: list[Allele] = []
    if x_count >= min_reads:
        alleles.append(_sex_allele(0, "X", x_count))
    if y_count >= min_reads:
        alleles.append(_sex_allele(1, "Y", y_count))

    if not alleles:
        return _no_data(system)

    call_rule = CallRule.HOMOZYGOUS if len(alleles) == 1 else CallRule.HETEROZYGOUS
    return MarkerResult(
        marker_name=system.name,
        system=system,
        alleles=alleles,
        alleles_called=alleles,
        call_rule=call_rule,
        tri_type=TriType.NONE,
        total_reads=total_reads,
    )


def _sex_allele(idx: int, letter: str, n_reads: int) -> Allele:
    return Allele(
        cluster_index=idx,
        consensus=letter,
        length_bp=0,
        n_reads_total=n_reads,
        n_reads_hp1=0,
        n_reads_hp2=0,
        n_reads_hp_none=n_reads,
        n_forward=0,
        n_reverse=0,
        mean_qual=0.0,
        ce=None,
        isfg=letter,
        bp_diff=0,
        is_deletion=False,
        status=AlleleStatus.ALLELE,
    )


def _no_data(system: System) -> MarkerResult:
    return MarkerResult(
        marker_name=system.name,
        system=system,
        alleles=[],
        alleles_called=[],
        call_rule=CallRule.NO_DATA,
        tri_type=TriType.NONE,
        total_reads=0,
    )


def _safe_pileup(
    bam: Path, chrom: str, start: int, end: int, min_mapq: int,
    *, reference_fasta: Path | None = None,
) -> list:
    try:
        return pileup_locus(bam, chrom, start, end, min_mapq=min_mapq,
                            reference_fasta=reference_fasta)
    except Exception:
        return []
