# FRONTStr

> **F**orensic **R**anked **O**utput for **N**anopore **T**andem **S**hort **T**andem **R**epeats

A forensic STR profiling toolkit for **long-read sequencing**, built ONT-first:
every default, calibration and threshold in it was measured on Oxford Nanopore
R10 data rather than carried over from capillary electrophoresis or short-read
practice.

That is the point of the project. A long read spans an entire STR locus and its
flanks at once, which makes it possible to report the **sequence** of each
allele rather than only its length: iso-alleles, microvariants and interrupted
repeat structures that CE cannot resolve. But a nanopore read's error profile
is nothing like a capillary trace's, and a caller that inherits CE-era
constants is confidently wrong in ways that are hard to see. Where FRONTStr had
inherited such constants, they have been replaced with values fitted to real
ONT data, and the measurements are in [`docs/`](docs/), including the ones that
changed the answer.

---

## Contents

1. [Status](#status)
2. [What it does](#what-it-does)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Quick start](#quick-start)
6. [Tutorial: from a public sample to a report](#tutorial-from-a-public-sample-to-a-report)
7. [Commands](#commands)
8. [Parameters](#parameters)
9. [Watching a run](#watching-a-run)
10. [Output formats](#output-formats)
11. [File formats](#file-formats)
12. [Cohort analysis](#cohort-analysis)
13. [Data requirements](#data-requirements)
14. [Limitations](#limitations)
15. [Documentation](#documentation)
16. [Development](#development)
17. [License and citation](#license-and-citation)

---

## Status

**Pre-alpha.** The pipeline works end to end from an aligned BAM and is
regression-tested against a real ONT sample. It has **not** been validated by a
forensic laboratory and is not fit for casework. The
[Limitations](#limitations) are specific and worth reading before you draw
conclusions from any output.

---

## What it does

**FRONTStr calls genotypes by itself.** Every allele, read count and sequence in
its output comes from its own pileup of the BAM.

| Stage | Module | What it does |
|---|---|---|
| Evidence | `frontstr.evidence` | Per-locus pileup → clustering by repeat-core length → POA consensus. **Integer per-allele read counts straight from the BAM.** |
| Interpretation | `frontstr.interp` | ISFG nomenclature, allele numbering, stutter, haplotype-aware suppression, genotype calling. |
| QC & audit | `frontstr.interp.qc`, `frontstr.audit` | Coverage, strand-bias and nomenclature flags; the sealed audit record. |

Coverage is counted, not divided out of a caller's `BPDIFFS` field, and the ISFG
bracket string is computed from the cluster consensus rather than from a VCF
`ALT`.

Four choices drive most of what the tool reports, each measured rather than
asserted:

- **POA consensus is required, not optional.** Without it the consensus is a
  single unpolished read, which on the test set manufactured four false
  microvariants in 202 called alleles.
- **Reads are binned by repeat-core length**, not window length. Panel windows
  are ~80% flank. Binning by repeat *unit count* was evaluated and rejected: it
  would merge TH01 9 with 9.3.
- **Stutter is log-linear in LUS** (R² 0.965), because the measured rate spans
  more than 10× across LUS 10–14.
- **Haplotype tags are evidence.** Two candidates on the same haplotype within
  2 bp cannot both be real alleles.

**→ [`docs/how_it_works.md`](docs/how_it_works.md) walks every stage with worked
examples, the derivation behind each threshold, and the full current state.**

---

## Requirements

### Platforms

| Platform | Status | Notes |
|---|---|---|
| **Linux** x86_64 | ✅ Fully supported | Binary wheels exist for every dependency. Simplest install. |
| **macOS** Intel & Apple Silicon | ✅ Supported | `pyspoa` has no macOS wheel and is compiled from source, so a C++ toolchain and CMake are needed. See below. |
| **Windows** | ⚠️ **Via WSL2 only** | Not a FRONTStr limitation. See below. |

#### Why Windows needs WSL2

FRONTStr reads BAM/CRAM through **pysam**, which wraps **htslib**. htslib is a
POSIX C library and pysam publishes **no Windows wheels and does not support
Windows**. The same is true of `cyvcf2` (VCF), `edlib` (alignment) and `pyspoa`
(POA consensus). Checked against PyPI on 2026-07-30: of the five compiled
dependencies, **zero** ship a `win_amd64` wheel.

This is upstream of FRONTStr and cannot be worked around by packaging. Every
htslib-based tool is in the same position, samtools included.

**WSL2 is the supported path and it is a full Linux environment**, so
performance and behaviour match native Linux:

```powershell
wsl --install -d Ubuntu
```

Then follow the [Linux instructions](#linux) inside the WSL2 shell. Your Windows
drives are mounted under `/mnt/c/…`, so BAMs on the Windows side are readable
without copying.

### Software

- **Python ≥ 3.11** (3.11, 3.12 and 3.13 are tested)
- **A POA backend.** Required, not optional. See [Installation](#installation).
- **samtools.** Not needed by FRONTStr itself, but needed to index a BAM and to
  cut region slices. Almost every workflow uses it.
- **minimap2.** Only if you are starting from FASTQ. FRONTStr takes aligned
  reads; the alignment step is not wired in.

---

## Installation

FRONTStr is installed **from a clone**, not from PyPI. The shipped panel and the
example files live in `examples/`, which is not part of the Python package, so a
wheel install would leave you without a panel to run against.

### Linux

```bash
git clone https://github.com/Tfronta/FRONTStr.git
cd FRONTStr
python3 -m venv .venv && source .venv/bin/activate
pip install -e '.[poa]'
```

### macOS

`pyspoa` publishes only Linux wheels, so pip builds it from source. That needs
CMake and the Apple command-line tools:

```bash
xcode-select --install          # once per machine
brew install cmake samtools
```

then the same three commands as Linux:

```bash
git clone https://github.com/Tfronta/FRONTStr.git
cd FRONTStr
python3 -m venv .venv && source .venv/bin/activate
pip install -e '.[poa]'
```

> **Apple Silicon:** use `pyspoa` (the `[poa]` extra above). The alternative
> backend `pyabpoa` hardcodes AVX2 x86 intrinsics and does not build on arm64.

### Windows

Install WSL2 (above), then follow the Linux instructions inside it.

### Verify the installation

```bash
frontstr doctor
```

This checks the environment only and needs no data. It must report a POA
backend. **A run without one is a different run**: the consensus degrades to a
single unpolished read and manufactures microvariants that are not in the
sample.

```
FRONTStr          0.1.0.dev0
Python            3.13.14
Platform          Linux x86_64
POA backend       poa_spoa
STRNaming         ready, 23 markers in the bundled slice cache
pysam             0.24.0
  ✓ installation looks complete
```

To also check that a BAM and the panel agree on contig names and coordinates:

```bash
frontstr doctor --bam sample.bam --panel examples/panels/codis_20_grch38.yaml
```

### Optional extras

| Extra | Install | Adds |
|---|---|---|
| `poa` | `pip install -e '.[poa]'` | **Required.** POA consensus via `pyspoa` |
| `poa-abpoa` | `pip install -e '.[poa-abpoa]'` | Alternative POA backend, Linux/Intel only |
| `parquet` | `pip install -e '.[parquet]'` | Parquet output for cohort tidy datasets |
| `pdf` | `pip install -e '.[pdf]'` | PDF rendering of the HTML report |
| `dev` | `pip install -e '.[dev,poa]'` | Test suite, linters, type checker |

---

## Quick start

A real ONT sample ships with the repository, so this runs on a fresh clone with
**nothing else downloaded**:

```bash
frontstr interpret --bam demodata/HG00113.demo.bam --panel examples/panels/codis_20_grch38.yaml
```

```
┏━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━┳━━━━━━━━━━━┓
┃ Marker   ┃ Call         ┃ Tri ┃ Alleles called ┃   AB ┃ QC        ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━╇━━━━━━━━━━━┩
│ D3S1358  │ homozygous   │ -   │ 14(31)         │    - │ -         │
│ vWA      │ heterozygous │ -   │ 14(17), 16(16) │ 0.52 │ -         │
│ FGA      │ heterozygous │ -   │ 24(15), 21(11) │ 0.58 │ -         │
```

Then, in order of how much they tell you:

```bash
# how it reached every call: reads, bins, clusters, haplotypes, naming
frontstr interpret --bam demodata/HG00113.demo.bam --panel examples/panels/codis_20_grch38.yaml --trace

# the full set of output files, including a self-contained HTML report
frontstr export --bam demodata/HG00113.demo.bam --panel examples/panels/codis_20_grch38.yaml \
  --out-dir out/ --formats profile,seqs,evidence,json,html,xlsx
```

`demodata/HG00113.demo.bam` is 17 MB: 1,077 ONT R10 reads from a public
1000 Genomes sample, covering all 25 panel markers. It reproduces the reference
profile exactly. See [`demodata/README.md`](demodata/README.md) for its
provenance and expected output.

To run on **your own** data, point `--bam` at any indexed BAM or CRAM aligned to
GRCh38. The [tutorial](#tutorial-from-a-public-sample-to-a-report) below shows
how to cut a panel-sized slice from a whole-genome BAM.

---

## Tutorial: from a public sample to a report

The [quick start](#quick-start) already ran on the bundled demo. This section
does the step the demo skips: **cutting a panel-sized slice out of a
whole-genome BAM**, which is what you will do with your own samples. It uses
**HG00113** from the public ONT open-data bucket, the same sample the demo was
made from, so you can check your result against a known answer.

### 1. Cut a panel-sized slice from the public BAM

The published BAMs are whole-genome and tens of gigabytes. Nothing needs to be
downloaded whole: samtools is pointed at the remote URL with its remote index
and asked only for the panel regions, which pulls the byte ranges those regions
occupy and nothing else.

```bash
BAM=https://s3.amazonaws.com/1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/HG00113-ONT-hg38-R10-LSK114-dorado090_sup_5mCG_5hmCG_v500.phased.bam
samtools view -b -M -L tests/data/ont_slices/codis_pm10kb.bed -X "$BAM" "$BAM.bai" -o HG00113.bam
samtools index HG00113.bam
```

About **31 MB** and **under a minute** on a normal connection. `-M` makes `-L` a
true multi-region filter; `-X` names the remote index so htslib never scans the
whole object.

### 2. Check the environment and the BAM against the panel

```bash
frontstr doctor --bam HG00113.bam --panel examples/panels/codis_20_grch38.yaml
```

### 3. Call the profile

```bash
frontstr interpret --bam HG00113.bam --panel examples/panels/codis_20_grch38.yaml
```

```
┏━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━┳━━━━━━━━━━━┓
┃ Marker   ┃ Call         ┃ Tri ┃ Alleles called ┃   AB ┃ QC        ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━╇━━━━━━━━━━━┩
│ D3S1358  │ homozygous   │ -   │ 14(31)         │    - │ -         │
│ vWA      │ heterozygous │ -   │ 14(17), 16(16) │ 0.52 │ -         │
│ FGA      │ heterozygous │ -   │ 24(15), 21(11) │ 0.58 │ -         │
```

### 4. Follow one locus back to the reads

```bash
frontstr interpret --bam HG00113.bam --panel examples/panels/codis_20_grch38.yaml --trace
```

Per locus this prints the read funnel with a named reason for every rejected
read, the repeat-core bins, each cluster with its consensus and haplotype
counts, how each allele number was derived, and why each discarded candidate was
discarded:

```
── vWA  chr12:5983877-5984149  (273 bp window, motif TCTA,TCTG, compound, minus strand)
  Reads fetched around the window         45
  Rejected (7)
      not a primary alignment             3
      MAPQ below threshold                0
      does not reach the left flank anchor  4
      ...
  Spanning the whole window               38   (total locus coverage)

  Step 1: grouped by length              7 bins, using repeat-core length
      56 bp core                          17 reads
      72 bp core                          16 reads
  Step 2: split by sequence              7 clusters
  Step 3: consensus per cluster          poa_spoa

  Candidates, strongest first
    * #0    17 reads  44.7%   261 bp  HP1 0 / HP2 17 / untagged 0  [block 4552958]
        name      CE14_GGAT[5]AGAT[3]GGAT[1]AGAT[2]AGAC[4]AGAT[1]AGAC[1]AGAT[2]
        number    14   via STRNaming
        verdict   allele
    ...
  Genotype                                14 (17 reads), 16 (16 reads)   [heterozygous]
  Allele balance                          0.52  (balanced)
```

To narrow it to one locus:

```bash
frontstr interpret --bam HG00113.bam --panel examples/panels/codis_20_grch38.yaml --trace 2>&1 >/dev/null \
  | awk '/^── vWA /{f=1} f&&/^── /&&!/vWA/{exit} f'
```

### 5. Write the full set of outputs

```bash
frontstr export --bam HG00113.bam --panel examples/panels/codis_20_grch38.yaml \
  --out-dir out/ --formats profile,seqs,evidence,json,html,xlsx \
  --operator "your name" --run-id RUN-001
```

Open `out/HG00113.html` in a browser: the forensic profile with ISFG
nomenclature, the sequence view, QC, per-locus detail, and the audit record.

### 6. Run several samples

Write a tab-separated manifest:

```
sample_id	bam	role
HG00113	/data/HG00113.bam	sample
HG00114	/data/HG00114.bam	sample
CTRL_POS	/data/control.bam	positive_ctrl
```

```bash
frontstr batch --manifest samples.tsv --panel examples/panels/codis_20_grch38.yaml \
  --out out/ --formats profile,evidence,seqs,json,html -j 4
```

This writes a per-sample directory, `batch_summary.csv`, a tidy cohort dataset,
and, for more than one sample, **`cohort.html`**: every sample at every marker,
one marker per tab, with a Download Excel button and a link from each sample to
its own report. See [Cohort analysis](#cohort-analysis).

---

## Commands

| Command | What it does |
|---|---|
| `interpret` | **The canonical command.** Pileup → cluster → ISFG → classify → call, printed as a table |
| `export` | The same run, written to files in the formats you name |
| `report` | The same run, written as a self-contained HTML report |
| `batch` | Many samples from a manifest, in parallel, plus the cohort view |
| `evidence` | Per-locus cluster dump. A debugging view of the evidence layer |
| `tidy` | Build a cohort-scale long table from run JSONs |
| `calibrate-stutter` | Fit a stutter model to your own data |
| `doctor` | Pre-flight: the environment, then BAM ↔ panel compatibility |
| `inspect` | Detect and validate an input file |
| `run` | **Not implemented.** The FASTQ → alignment → report path is a stub |

`frontstr <command> --help` documents every flag.

---

## Parameters

Ten parameters decide every call. They are printed at the start of every run,
because the values that decide a genotype are usually the ones nobody typed.

| Parameter | Default | Provenance | What it decides |
|---|---|---|---|
| `--min-mapq` | 20 | chosen | Reads below this MAPQ are dropped. Also what excludes X/Y paralogue mismappings at the sex markers |
| `--flank-anchor` | 20 | chosen | bp of cleanly aligned flank required each side. Keeps partially spanning reads, whose repeat tract is truncated, out of the pileup |
| `--identity` | 0.97 | chosen ⚠️ | Pairwise identity to join a cluster. **Known miscalibrated**, see the docs |
| `--len-tolerance` | 0 | **derived** | Merges adjacent length bins. **Must stay 0**: any tolerance merges TH01 9 with 9.3 |
| `--analytical-thresh` | 0.02 | chosen | Below this fraction of locus coverage, a cluster is noise and is not shown |
| `--calling-thresh` | 0.10 | chosen | Below this fraction, a cluster is an artefact and is shown but not called |
| `--min-phr` | 0.40 | convention | Minor allele as a fraction of the major before a heterozygote is called on read counts. Inherited from CE, and overridden by haplotype evidence |
| `--min-reads-third` | 5 | **derived** | Absolute read floor a 3rd candidate must clear before it can make a locus triallelic or raise a mixture flag |
| `--low-coverage-reads` | 20 | **derived** | A called locus below this many supporting reads is flagged |
| `--balanced-ab-max` | 0.65 | chosen | Largest allele balance a heterozygote may have and still count as balanced |

**Provenance is part of the value.** `derived` means computed from measured data;
overriding one of those raises `NON_DEFAULT_THRESHOLD` on **every marker**, so
the profile still says six months later which numbers were changed and which had
measured backing. `chosen` means defensible but picked. `convention` means
inherited from forensic practice rather than from this data.

Each default's derivation is documented where the value lives, next to the
constant it justifies: `frontstr/params.py` for the table, and
`frontstr/interp/qc.py` for the coverage floor and the balance band. That is
deliberate, so a number and its reasoning cannot drift apart.

---

## Watching a run

Two flags make a run inspectable, and both write to **stderr**, so redirecting
the results table still separates cleanly.

| Flag | Gives you |
|---|---|
| `--log` | The configuration, then one line per marker: call rule, coverage, cluster count, how each allele number was derived |
| `--trace` | The full per-locus narrative shown in the tutorial above |
| `--trace-out FILE` | The narrative to a file instead of the terminal |

`frontstr batch` writes `<out>/<sample>/<sample>.trace.txt` **by default**, one
file per sample. That is deliberate: the trace is the record that lets someone
ask *where* a call went wrong rather than only whether it did, and it costs
nothing measurable (18.5 s against 19.0 s over five samples, ~160 kB each).
`--no-trace` exists for a constrained output directory; it is not a path anyone
should have to remember. With `-j 1` the narrative also streams to the terminal
as it happens.

---

## Output formats

| `--formats` | Output |
|---|---|
| `profile` | CSV, one row per marker. The genotype table |
| `evidence` | CSV, one row per cluster, including the uncalled ones |
| `seqs` | CSV, one row per called allele with ISFG and consensus |
| `json` / `json-compact` | The canonical record: everything, plus the audit block |
| `html` | Self-contained offline report |
| `xlsx` | Five-sheet review workbook (Profile, Sequences, Evidence, QC, Audit) |
| `vcf` | Native sequence-resolved VCF. **Needs `--reference`** |

The VCF is native, not an annotation of another caller's output. `ALT` carries
the allele's sequence so iso-alleles stay distinct; the repeat count rides in
`FORMAT/MC` and per-allele depth in `FORMAT/AD`. It is bgzip/tabix indexable and
bcftools-queryable, which is the point: it exists so FRONTStr can be benchmarked
against other callers.

Alongside the exports every run writes a `frontstr.log.jsonl` process log, and
the canonical JSON carries a sealed audit record: tool version, POA backend,
stutter model, every threshold that moved a call, input hashes, the flag census
and a SHA-256 over the record itself.

### Reading the profile table

| Column | Meaning |
|---|---|
| `Allele 1` / `Allele 2` | The called allele numbers |
| `AD` | Allelic depth: reads supporting that allele |
| `AB` | Allele balance: the strongest called allele over the called pair. 0.50 is even, 1.0 is everything on one allele. Only defined for a heterozygote |
| `QC` | Abbreviated flag codes, expanded in the legend under the table |
| `ISFG` | The bracketed repeat structure, per allele |

**There is no locus total, on purpose.** Depth is reported per allele and only
per allele, because that is the evidence standing behind each allele
separately, and it is what clustering the reads by sequence buys you. A caller
that only sizes the repeat cannot produce it. The spanning read count still
exists, since it is the denominator every fraction threshold is measured
against, and it is reported in the VCF as `DP`.

---

## File formats

### Region file (BED)

`--bed` points FRONTStr at your own regions without curating a panel YAML. The
file is a plain, headerless, tab-separated BED.

**NOTE: the table header below is descriptive. The BED file must not have one.**

| CHROM | START | END | MOTIF | NAME |
|---|---|---|---|---|
| chr1 | 230769515 | 230769783 | TAGA,CAGA | D1S1656 |
| chr10 | 129294140 | 129294395 | GGAA | D10S1248 |
| chr11 | 2170987 | 2171215 | AATG | TH01 |
| chr12 | 5983876 | 5984149 | TCTA,TCTG | vWA |
| chr2 | 1489550 | 1489684 | AATG | TPOX |

| Column | Required | Contents |
|---|---|---|
| 1 CHROM | yes | Contig name, matching the BAM header exactly |
| 2 START | yes | Window start. 0-based half-open by default; `--bed-coords panel1` reads it as 1-based inclusive |
| 3 END | yes | Window end |
| 4 MOTIF | **yes** | Repeat motif on the canonical strand. Comma-separated for compound markers |
| 5 NAME | no | Marker name. Defaults to `chrom_start_end` |

**The motif column is required, and that is not a parsing convenience.** Without
it reads cannot be binned by repeat-core length, which is the mechanism that
stops ONT flank indel errors splitting one allele into a dozen clusters. Binning
on raw window length instead took TH01 from 2 bins to 12. FRONTStr refuses the
file rather than calling badly from it.

A ready-made file for the shipped panel is in the repository:

```bash
frontstr interpret --bam demodata/HG00113.demo.bam --bed examples/panels/codis_20_grch38.bed
```

**A BED cannot carry calibration.** `corr_value`, `reference_ce`, strand and the
per-marker triallelic settings live in the panel YAML and have no BED column.
Markers STRNaming has no reporting range for therefore get an uncalibrated
repeat count instead of a kit allele, the run says so per locus, and every call
at them raises `CE_NOMENCLATURE_OFFSET`. Use the YAML panel for the CODIS loci
and a BED for exploratory regions.

> `demodata/demo_regions.bed` is **not** usable with `--bed`. It has no motif
> column because its only job is to tell `samtools` which regions to cut.

### VCF output

`--formats vcf` writes a native, sequence-resolved VCF. `ALT` carries each
allele's full sequence, so iso-alleles stay distinct instead of collapsing to a
repeat count. It needs `--reference`.

**INFO fields** describe the locus.

| Field | Number | Type | Description |
|---|---|---|---|
| `END` | 1 | Integer | End of the panel window |
| `MARKER` | 1 | String | Forensic marker name |
| `MOTIF` | 1 | String | Repeat motif(s), comma-separated, on the canonical strand |
| `PERIOD` | 1 | Integer | Repeat period in bp. `-1` for compound markers whose allele number comes from the bracket count |
| `STRAND` | 1 | String | Strand of the canonical ISFG motif relative to the reference |
| `NOTE` | . | String | Informational flags that do not filter the call |

**FORMAT fields** describe the genotype. Fields with `Number=R` carry one value
per allele in REF-then-ALT order.

| Field | Number | Type | Description |
|---|---|---|---|
| `GT` | 1 | String | Genotype. Triallelic loci are emitted as three alleles, e.g. `1/2/3` |
| `DP` | 1 | Integer | Reads spanning the whole window at this locus |
| `AD` | R | Integer | Allelic depth: reads supporting each allele |
| `AL` | R | Integer | Allele length in bp |
| `MC` | R | String | Canonical allele number, the CE or repeat count |
| `MM` | R | String | How each `MC` was derived: `strnaming`, `period_ce`, `bracket_count`, … |
| `ISFG` | R | String | ISFG bracketed nomenclature |
| `ISO` | R | String | Iso-allele suffix from the curated catalog, when one applies |
| `HP` | R | String | Read partition as `hp1｜hp2｜untagged` |
| `SB` | R | String | Strand partition as `forward｜reverse` |
| `CM` | R | String | Consensus method: `poa_spoa`, `poa_abpoa` or `single` |

**FILTER** carries the QC flags. `PASS` means no warning or error flag was
raised; every other code is one of the QC conditions, so a VCF can be filtered
on exactly what the report shows.

> **`DP` appears only in the VCF.** There it is every read spanning the window,
> which is what the specification means by depth and what `bcftools` users
> expect. The profile tables report no locus total at all, only per-allele `AD`,
> so the two can never be confused for one another.

## Cohort analysis

For more than one sample, `frontstr batch` writes **`cohort.html`** beside the
per-sample reports: one block per marker, one row per sample, with both alleles,
their depths, allele balance, QC and both ISFG strings.

- **One marker per tab.** 108 samples across 25 markers stacked would be 2,700
  rows of continuous scroll. Each tab carries the marker's call rate or flagged
  count, so the tab strip itself says which locus to open.
- **"Where to look first"** names the five worst markers and the five worst
  samples. A locus that fails across the whole cohort is a statement about the
  panel or the assay, not about the samples.
- **Download Excel** writes a real workbook: `Genotypes` in long form (one row
  per marker × sample, filterable and pivotable), plus `By marker` and
  `By sample` margins.
- **Each sample name links to its own report**, and that report links back.

For analysis rather than review, `frontstr tidy` flattens run JSONs into one long
table, one row per **sample × marker × allele**, as CSV and Parquet:

```bash
frontstr tidy --from-dir out/ -o analysis/
```

`batch` produces it automatically when `json` is among its formats. Because it is
built from the canonical JSONs, a dataset can be rebuilt at any time without
re-running, and runs from different batches combined. Markers that dropped out
appear as `called = false` rows rather than vanishing, and each row carries the
panel version, POA backend and stutter model that produced it, so a cohort
collected across two calibrations can still be analysed honestly.

---

## Data requirements

### Input reads

- **Aligned, coordinate-sorted, indexed BAM or CRAM.** CRAM additionally needs
  `--reference`.
- **Aligned to GRCh38** for the shipped panel. `minimap2 -ax map-ont` against an
  hg38 reference is what the panel's coordinates and calibration assume. The
  reference profiles were produced against `hg38.no_alt`.
- **ONT R10 chemistry, Dorado basecalling.** The stutter model and every
  `corr_value` in the panel were calibrated on R10 + Dorado. R9 or guppy data
  will run but is being scored against a model that does not describe it.
- **Phased reads are optional but valuable.** `HP`/`PS` tags (from longphase,
  whatshap or similar) enable two calling rules: haplotype-aware phantom
  suppression and the peak-ratio rescue. Both are no-ops on an unphased BAM.

### The panel

`examples/panels/codis_20_grch38.yaml` holds the CODIS 20 core loci plus AMEL and
four sex-validation STRs, 25 systems in total, with GRCh38 coordinates, motifs
and calibrated `corr_value` per marker.

To use your own regions instead, `--bed` takes a plain BED. It cannot carry
calibration, so markers STRNaming has no reporting range for get an uncalibrated
repeat count rather than a kit allele, and the run says so per locus.

### Example data

One sample **is** versioned: [`demodata/HG00113.demo.bam`](demodata/README.md),
17 MB, so a fresh clone can run the tool without downloading anything. It is a
±1 kb cut of a public 1000 Genomes ONT BAM with the unused tags stripped, and it
reproduces the reference profile exactly.

The wider ±10 kb slices the full regression suite uses (five samples, ~200 MB)
are **not** versioned. Tests needing them skip cleanly; the HG00113 assertions
run against the demo instead. Build your own with the tutorial's step 1, which
is how they were made.

---

## Limitations

Part of the documentation, not a disclaimer.

- **`frontstr run` is a stub.** The FASTQ → alignment → report path is not
  wired; `ingest.align` raises. Align externally and start from a BAM.
- **Not laboratory-validated.** No forensic partner has signed off. No mixture
  series, no dropout study, no NIST control. The reference profile it is tested
  against comes from another caller on matched Illumina data. That is
  caller-vs-caller concordance, not validation against a reference method.
- **The stutter model is PCR-free.** Fitted on WGS, so it has no PCR slippage
  component and will under-predict stutter on an amplicon panel. The model
  records its own protocol so this cannot be lost.
- **Two markers still have no ISFG designation.** DYS393 is absent from
  STRNaming's reported-range table and AMEL is not a tandem repeat, so both keep
  the legacy length/bracket arithmetic. `Allele.number_method` records which path
  produced each number.
- **The analytical and calling thresholds (0.02 / 0.10) are not data-derived.**
  The coverage floor and the stutter model are; these two are chosen defaults,
  and that is a gap.
- **`--identity` is known to be miscalibrated.** It was set as if comparing a
  read to a consensus, but the comparison is raw read to raw read, where two ONT
  reads of one allele diverge 2–4%.
- **No CODIS `.cmf`, MIDST or PDF export.**
- **Illumina data will not work** against the shipped panel: its windows are
  ±100 bp and the pileup needs reads that span them fully. A deliberate ONT-first
  choice, not an oversight.

### Planned

- Validation across the published 1000 Genomes ONT cohort, against Illumina
  genotypes for the same samples
- FASTQ input, with the alignment step wired into `frontstr run`
- A stutter model fitted on amplicon data, alongside the PCR-free one
- A per-locus aligned-read view, to inspect the pileup a call was made from
- CODIS `.cmf` export

---

## Documentation

| Document | Contents |
|---|---|
| [`docs/how_it_works.md`](docs/how_it_works.md) | Every stage step by step, with examples and measurements; full current state |
| [`docs/stutter_calibration.md`](docs/stutter_calibration.md) | How the stutter model was measured, and its caveats |
| [`docs/architecture.md`](docs/architecture.md) | Module layout |
| [`demodata/README.md`](demodata/README.md) | The bundled demo sample: provenance, how it was cut, and its expected profile |
| [`benchmark/README.md`](benchmark/README.md) | The development benchmark against other callers. **Not part of running FRONTStr** |

---

## Development

```bash
pip install -e '.[dev,poa]'
pytest                      # unit suite
pytest -m integration       # end-to-end regression, needs local ONT slices
ruff check frontstr tests && ruff format --check frontstr tests && mypy frontstr
```

The integration tests run the whole pipeline against a real ONT R10 BAM and
assert the genotype of every marker. **They run on a fresh clone**, against the
bundled [demo sample](demodata/README.md), so the end-to-end assertions are
exercised in CI rather than only on one laptop. The tests that need the other
four samples skip cleanly until those ~200 MB slices are rebuilt. See
[Data requirements](#data-requirements).

---

## License and citation

MIT. See [`LICENSE`](LICENSE).

If you use FRONTStr in research, please cite:

> [TBD, pending preprint]

### Prior art

FRONTStr is an independent implementation. It uses and draws on:

- **LongTR** [10.1186/s13059-024-03319-2](https://doi.org/10.1186/s13059-024-03319-2)
- **ISFG sequence nomenclature.** The bracketed-repeat convention FRONTStr emits
  follows the ISFG recommendations and the counting rules described by Gettings
  et al. (2024) [10.1016/j.fsigen.2023.102946](https://doi.org/10.1016/j.fsigen.2023.102946)
- **toaSTR.** Ganschow S, Silvery J, Kalinowski J, Tiemann C (2018)
  [10.1016/j.fsigen.2018.07.006](https://doi.org/10.1016/j.fsigen.2018.07.006)
- **STRspy.** Hall, Kesharwani et al.
  [10.3390/ijms27041889](https://doi.org/10.3390/ijms27041889)
