# Cohort validation harness

Runs FRONTStr over the 1000 Genomes ONT cohort and scores it against the
external truth genotypes. Dev tooling — `pyproject.toml` ships only `frontstr`,
so nothing here is installed with the package and nothing in the caller imports
it.

## Why it is not a `frontstr` command

The truth workbook is not in this repository and its layout is bespoke: merged
header rows, Spanish field names, three technologies side by side, marker blocks
of unequal width. Baking that into the shipped CLI would make a private
spreadsheet part of the package's contract. The seam is deliberate —
`truth.py` knows about the workbook, and nothing downstream of it does.

## The three steps

```bash
python -m validation.cohort fetch --out cohort/slices \
    --workbook ~/Desktop/1000GEN-ONT-Merged-Compar.xlsx

frontstr batch --manifest cohort/slices/manifest.tsv \
    -p examples/panels/codis_20_grch38.yaml -o cohort/run \
    --formats profile,json -j 4

python -m validation.compare --summary cohort/run/batch_summary.csv \
    --workbook ~/Desktop/1000GEN-ONT-Merged-Compar.xlsx --out cohort
```

`python -m validation.cohort plan` reports what `fetch` would download without
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
