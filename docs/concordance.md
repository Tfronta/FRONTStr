# Concordance against Illumina, 108 ONT samples

A measurement, not a validation. It records how often FRONTStr's genotype
matches an Illumina genotype for the same individual, on a public cohort. It
does not establish that either is correct, and no laboratory has validated
either method here. Read it alongside [Limitations](../README.md#limitations).

## The cohort

108 samples from the 1000 Genomes ONT release, every one R10 chemistry
basecalled with Dorado. R9 and guppy calls are rejected: the stutter model and
the panel's `corr_value` table were measured on R10 + Dorado, so an R9 sample
would be scored against a model that does not describe it.

Panel windows only, sliced from the public bucket with `samtools` against the
remote index, roughly 44 MB a sample. The harness that does this is
[`benchmark/`](../benchmark/README.md) and is not part of the caller.

## Result

FRONTStr 0.1.0.dev0, panel `codis_20_grch38`, default parameters.

| | |
|---|---|
| Loci called | 2160 |
| Loci with an Illumina genotype to compare against | 2045 |
| Match | 1997 |
| Mismatch | 48 |
| **No call** | **0** |
| Concordance | **97.7 %** |

Per marker:

| Marker | Scored | Mismatch | Concordance |
|---|---|---|---|
| D10S1248 | 108 | 0 | 100.0 % |
| D2S1338 | 108 | 0 | 100.0 % |
| TPOX | 108 | 0 | 100.0 % |
| CSF1PO | 108 | 1 | 99.1 % |
| D13S317 | 107 | 1 | 99.1 % |
| D16S539 | 108 | 1 | 99.1 % |
| D22S1045 | 108 | 1 | 99.1 % |
| TH01 | 108 | 1 | 99.1 % |
| D12S391 | 107 | 2 | 98.1 % |
| D3S1358 | 108 | 2 | 98.1 % |
| D5S818 | 108 | 2 | 98.1 % |
| D7S820 | 106 | 2 | 98.1 % |
| D8S1179 | 108 | 2 | 98.1 % |
| vWA | 108 | 2 | 98.1 % |
| D19S433 | 108 | 3 | 97.2 % |
| D2S441 | 108 | 4 | 96.3 % |
| D1S1656 | 107 | 5 | 95.3 % |
| FGA | 107 | 8 | 92.5 % |
| D18S51 | 107 | 11 | 89.7 % |
| D21S11 | 0 | — | no Illumina genotype in the comparison set |

**115 loci are excluded for want of a comparison, not skipped.** 108 of them are
D21S11, which the comparison set leaves blank for every sample. AMEL and the
four sex markers (DYS391, DYS393, DXS7132, DXS8378) are absent from it
altogether and are not scored at all, so nothing in this table speaks to them.

## What the 48 mismatches are

Classified by mechanism from the per-locus evidence, not by which caller was
assumed right. The counts are from the run before the identical-consensus merge
(50 mismatches); the merge changed four loci, listed under it.

**Three sequence callers against one length caller — 11 loci.** FRONTStr,
longTR and STRspy agree with each other and differ from Illumina. Both other
callers read the same ONT reads, so this is not independent confirmation of
FRONTStr, but it does locate the disagreement between technologies rather than
in any one caller.

**Microvariant designation — 9 loci.** The integer repeat count agrees and the
fractional suffix does not. Four are D2S441, all with the same shape: FRONTStr
reports `.1`, Illumina `.2`, and STRspy agrees with FRONTStr in all four. A
question about naming conventions, answerable against the ISFG recommendation
rather than by counting agreements.

**Allele dropout — 14 loci, 0.68 % of those scored.** FRONTStr calls homozygous
and both other callers find the second allele in the same reads. This is the
one class that is a defect in FRONTStr. Spread over D1S1656 (4), D18S51 (3),
D5S818 (2), FGA (2), D8S1179 (2) and TH01 (1).

Two mechanisms, from the per-locus evidence:

- *Rejected by the heterozygote ratio.* The second allele is a called
  candidate and does not reach `min_phr_for_het`. In five of these the minor
  allele carries exactly two haplotype-tagged reads and
  `haplotype.DEFAULT_MIN_TAGGED_READS` is three, so the phasing rescue cannot
  fire; in two more the major cluster's haplotype purity is 0.75 or 0.78
  against a `DEFAULT_HP_PURITY` of 0.80.
- *Split below the ratio.* One allele occupies more than one cluster, and no
  fragment reaches the ratio on its own.

**Genotypes carrying more than two alleles — 3 loci.** Two of them repeat a
number (`10/12/12`, `21/22/22`), which is the split-cluster mechanism with both
fragments strong instead of both weak.

**Remaining — 13 loci.** No shared mechanism.

### Change from the identical-consensus merge

Merging clusters whose polished consensus is byte-identical
(`evidence.cluster._merge_identical_consensus`) moved four loci and took the
total from 50 to 48.

| Sample | Marker | Before | After | Illumina |
|---|---|---|---|---|
| HG01816 | FGA | `20/20` | `20/23` | `20/23` |
| HG02187 | TPOX | `10/12/12` | `10/12` | `10/12` |
| HG04161 | FGA | `26/26` | `24.2/26` | `24.2/26` |
| HG02555 | D18S51 | `15/15` | `14.2/15` | `15/15` |
| HG03667 | FGA | `21/22/22` | `21/21.2/22` | `21/22` |

The HG02555 change is worth reading in full, because the earlier agreement was
not evidence of a correct call. The `14.2` differs from the `15` by one
substitution and one `AG` dinucleotide inside a run of `AG` adjacent to a long
`GAAA` tract, arrives on 1 forward and 7 reverse reads at a locus that is 12
and 11 overall, and has no haplotype supporting it. Before the merge its reads
sat in three clusters that each failed the heterozygote ratio, so the locus was
called homozygous for a reason unrelated to any of that. `STRAND_BIAS` says
nothing because its binomial has no power at eight reads. The merge did not
create this candidate; it assembled one that was already there.

## Reproducing it

Both inputs are outside the repository: a listing of the bucket's BAM paths,
and the workbook holding the Illumina, longTR and STRspy genotypes. See
[`benchmark/README.md`](../benchmark/README.md).

```bash
python -m benchmark.cohort fetch --out cohort/slices --workbook <workbook.xlsx>
frontstr batch --manifest cohort/slices/manifest.tsv \
    -p examples/panels/codis_20_grch38.yaml -o cohort/run \
    --formats profile,json -j 4
python -m benchmark.compare --summary cohort/run/batch_summary.csv \
    --workbook <workbook.xlsx> --out cohort
```

About two hours to fetch and eight minutes to call. The slices are kept, so
re-scoring after a change to the caller costs the eight minutes.

## How this number is and is not used

It measures. It is not an objective to raise.

Every threshold in [`params.py`](../frontstr/params.py) carries its provenance,
and `derived` means a value computed from measured data: `low_coverage_reads`
comes from a binomial dropout calculation, not from a concordance score. A
threshold moved because concordance improved would be fitted to this comparison
set, and the number would stop measuring anything.

The test that separates the two: **would the change stand if concordance had
fallen?** The identical-consensus merge would, because two clusters carrying
byte-identical sequence are one allele by definition of allele, and the cohort
found the case rather than justifying it. Loosening the phasing rescue's gates
would not, which is why they are unchanged.
