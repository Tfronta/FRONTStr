# Development benchmark

**This is not part of FRONTStr. Nothing here runs when you call a sample.**

FRONTStr is for anyone with their own ONT data, where there is no Illumina
profile, no truth table and no second caller to compare against. That is the
normal case, and the caller is built for it: a genotype has to be defensible
from the evidence behind it, not from agreeing with somebody else's answer.
`frontstr interpret --trace` is where that happens — the reads, the bins, the
clusters, the haplotypes and the naming, per locus.

This directory is the opposite situation, and it exists only for us. A handful
of public 1000 Genomes samples have Illumina and other-caller genotypes
published; scoring against them is how we answer "is the caller right" while
developing it. That is **benchmarking**, and it is kept under its own name so it
can never read as a step in calling a sample.

If you are running FRONTStr on your own samples, you want
[`docs/how_it_works.md`](../docs/how_it_works.md) and `--trace`. Not this.

## Why it is not a `frontstr` command

Two reasons, and the second is the important one:

1. The truth workbook is not in this repository and its layout is bespoke —
   merged header rows, Spanish field names, three technologies side by side,
   marker blocks of unequal width. Baking that into the shipped CLI would make a
   private spreadsheet part of the package's contract.
2. A caller that offers "compare against longTR" as a subcommand teaches the
   wrong thing: that a call is trustworthy when a second caller agrees. Most
   users have no second caller. The evidence view has to carry the weight.

## The three steps

```bash
python -m benchmark.cohort fetch --out cohort/slices \
    --workbook ~/Desktop/1000GEN-ONT-Merged-Compar.xlsx

frontstr batch --manifest cohort/slices/manifest.tsv \
    -p examples/panels/codis_20_grch38.yaml -o cohort/run \
    --formats profile,json -j 4

python -m benchmark.compare --summary cohort/run/batch_summary.csv \
    --workbook ~/Desktop/1000GEN-ONT-Merged-Compar.xlsx --out cohort
```

`python -m benchmark.cohort plan` reports what `fetch` would download without
downloading anything.

### Measured cost, on one sample

| | |
|---|---|
| Slice size | ~35 MB (panel windows ±10 kb, out of a whole-genome BAM) |
| Fetch time | ~63 s per sample, network-bound |
| Call time | ~4 s per sample; 5 samples in 19 s at `-j 4` |
| **108 samples** | **~2 h to fetch, ~4 GB on disk, ~8 min to call** |

`fetch` is resumable — a sample whose slice and index already exist is skipped,
and a download that dies leaves a `.partial` rather than a truncated BAM that
resume would accept. Keeping the slices is the point: re-scoring the cohort
after a caller change then costs the 8 minutes, not the 2 hours.

## What counts as validation

Illumina is the only external anchor. **longTR and STRspy are ONT callers on the
same reads**, so FRONTStr agreeing with them is caller-vs-caller, not
independent confirmation. Both are still printed on every discordance, because
"FRONTStr contradicts Illumina but matches both ONT callers" is a different
finding from "FRONTStr stands alone".

Markers with no external truth are excluded from the concordance percentage and
labelled in the table rather than silently counted:

- **D21S11** — `NA` in the Illumina column for every sample in the workbook.
- **The sex markers** (DYS391, DYS393, DXS7132, DXS8378) and AMEL are absent
  from this sheet entirely, so they are not scored at all.

Samples that were called but have no row in the workbook are reported by name,
not dropped — three of the five committed test slices have no truth, and a
silent join made a five-sample batch look like a two-sample validation.

## Pilot results

Three samples with truth (HG00112 fetched, HG00113 and HG00154 from the
committed slices): 57 loci scored, 54 match, **94.7 %**. The three discordances
are all HG00154 and all worth looking at rather than explaining away:

| Marker | FRONTStr | Illumina | longTR | STRspy |
|---|---|---|---|---|
| D18S51 | 14/14 | 13/14 | 13/14 | 14/14 |
| D5S818 | 13/13 | 12/13 | — | 12/13 |
| FGA | 19/24 | 19/22 | 19/24 | 19/24 |

The first two are FRONTStr calling homozygous where Illumina calls
heterozygous — allele dropout or the minor allele being absorbed as stutter. At
FGA both ONT callers agree with FRONTStr against Illumina, which is a
technology-level disagreement, not a FRONTStr defect.
