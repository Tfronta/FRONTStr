# Demo data

One real ONT sample, small enough to version, so that a fresh clone can run
FRONTStr end to end **without downloading anything**.

```bash
frontstr interpret --bam demodata/HG00113.demo.bam --panel examples/panels/codis_20_grch38.yaml
```

| File | Size | What it is |
|---|---|---|
| `HG00113.demo.bam` | 17 MB | 1,077 ONT reads covering the 25 panel markers |
| `HG00113.demo.bam.bai` | 745 KB | Index |
| `demo_regions.bed` | 794 B | The 26 regions it was cut from, panel windows ±1 kb |

## Where it comes from

**HG00113**, a British (GBR) individual from the **1000 Genomes Project**, as
released by the **1000 Genomes ONT Sequencing Consortium** in its public
open-data bucket:

```
s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/
  HG00113-ONT-hg38-R10-LSK114-dorado090_sup_5mCG_5hmCG_v500.phased.bam
```

Public, consented, unrestricted-use data. 1000 Genomes samples were collected
with broad consent for open release precisely so they can be redistributed in
cases like this one; nothing here needs a data-access agreement.

| Property | Value |
|---|---|
| Chemistry | ONT R10.4.1, LSK114 |
| Basecaller | Dorado 0.9.0, `sup` model |
| Reference | GRCh38, `hg38.no_alt`, minimap2 2.28 |
| Phasing | longphase 1.7.3 over Clair3 variants — `HP`/`PS` tags present |
| Sex | Male, which is why AMEL calls X,Y and DYS391/DYS393 are callable |
| Mean read length | ~14 kb |

## How it was made

Two things were done to the published BAM, and **neither changes a genotype**:

1. **Cut to the panel windows ±1 kb** instead of the ±10 kb the test slices use.
   ONT reads average 14 kb, so every read that spans a marker window is still
   here; only reads that merely touched the wider margin are gone. 3,504 reads →
   1,077.
2. **Dropped the tags FRONTStr does not read.** The caller uses only `HP` and
   `PS`; the methylation tags (`MM`/`ML`) and the aligner's own tags (`SA`,
   `NM`, `ms`, `AS`, …) are removed. Base qualities and sequences are untouched.

Reproduce it exactly:

```bash
BAM=https://s3.amazonaws.com/1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/HG00113-ONT-hg38-R10-LSK114-dorado090_sup_5mCG_5hmCG_v500.phased.bam

samtools view -b -M -L demodata/demo_regions.bed -X "$BAM" "$BAM.bai" \
  --output-fmt-option level=9 \
  -x MM,ML,MN,SA,ms,AS,nn,tp,cm,s1,s2,de,rl,PQ,NM \
  -o demodata/HG00113.demo.bam
samtools index demodata/HG00113.demo.bam
```

About 50 seconds and 17 MB. `demo_regions.bed` is derived from the panel, so it
can be regenerated if the panel windows ever move.

**Verified:** this file produces byte-for-byte the same 25 genotypes, allele
numbers and per-allele depths as the full ±10 kb slice. The size reduction costs
nothing.

## Expected output

The full profile, which `tests/test_regression_hg00113.py` asserts marker by
marker:

| Marker | Alleles | Depth | Marker | Alleles | Depth |
|---|---|---|---|---|---|
| D3S1358 | 14 | 31 | D19S433 | 14, 12 | 39 |
| vWA | 14, 16 | 33 | D10S1248 | 13, 15 | 30 |
| FGA | 24, 21 | 26 | D1S1656 | 13, 11 | 25 |
| D8S1179 | 10, 13 | 25 | D2S441 | 13, 10 | 32 |
| D21S11 | 31, 29 | 26 | D12S391 | 21, 20 | 33 |
| D18S51 | 14, 15 | 21 | D22S1045 | 16 | 22 |
| D5S818 | 11, 13 | 31 | AMEL | X, Y | 26 |
| D13S317 | 11, 9 | 23 | DYS391 | 10 | 14 |
| D7S820 | 8, 10 | 27 | DYS393 | 14 | 17 |
| D16S539 | 11, 8 | 33 | DXS7132 | 14 | 9 |
| TH01 | 9.3, 7 | 17 | DXS8378 | 10 | 10 |
| TPOX | 9, 11 | 35 | | | |
| CSF1PO | 10, 12 | 28 | | | |

Five markers raise QC flags, all of them `low_coverage` or phasing notes, and
all of them expected at this depth: TH01, DYS391, DYS393, DXS7132, DXS8378.
A quarter of loci flagged at ~30× ONT is the designed rate, not a symptom —
see [`docs/parameters.md`](../docs/parameters.md#low_coverage_reads--20--derived).

## What this sample is not

- **Not a validation set.** One sample, one individual, one run. It exists so
  the tool can be tried, not so it can be trusted.
- **Not a mixture, and not degraded.** Single-source, high-quality WGS. Nothing
  here exercises mixture detection or dropout.
- **Not amplicon.** PCR-free WGS, which is what the shipped stutter model was
  fitted on. An amplicon panel will stutter more than this model predicts.
