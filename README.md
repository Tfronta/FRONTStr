# FRONTStr

> **Forensic Ranked Output for Nanopore Tandem Short Tandem Repeats**

FRONTStr is a forensic STR profiling toolkit for **long-read sequencing data**,
built ONT-first: every default, every calibration and every threshold in it was
measured on Oxford Nanopore R10 data rather than carried over from capillary
electrophoresis or short-read practice.

That distinction is the point of the project. Long reads span an entire STR
locus plus its flanks in a single read, which makes it possible to report the
**sequence** of each allele rather than only its length — iso-alleles,
microvariants and interrupted repeat structures that CE cannot resolve. But the
error profile of a nanopore read is nothing like that of a capillary trace, and
a caller that inherits CE-era constants will be confidently wrong in ways that
are hard to see. Where FRONTStr had inherited such constants, they have been
replaced with values fitted to real ONT data — and the measurements are in
[`docs/`](docs/), including the ones that changed the answer.

## Status

**Pre-alpha.** The interpretation pipeline works end-to-end from an aligned BAM
and is regression-tested against a real ONT sample. It has **not** been
validated by a forensic laboratory, and it is not fit for casework. See
[Limitations](#limitations) — they are specific, and worth reading before you
draw conclusions from any output.

## What it does

Three layers, each with a strict responsibility
(see [`docs/architecture.md`](docs/architecture.md)):

| Layer | Module | Responsibility |
|---|---|---|
| 1 — Caller | `frontstr.caller` | Runs LongTR, parses its VCF. Used for cross-checking, never as the source of coverage. |
| 2 — Evidence | `frontstr.evidence` | Sequence-level pileup, clustering, POA consensus. **Integer per-allele read counts straight from the BAM.** |
| 3 — Interpretation | `frontstr.interp` | ISFG nomenclature, stutter, allele calling, QC flags. |

The forensic value sits in Layer 2: coverage is counted, not divided out of a
caller's `BPDIFFS` field, and the ISFG bracket string is computed from the
cluster's consensus rather than from a VCF `ALT`.

### Design decisions worth knowing about

These are the choices that most affect what the tool reports. Each is measured
and documented rather than asserted.

- **A POA consensus is required, not optional.** Without one the consensus
  degrades to a single unpolished read. Measured on five ONT R10 samples, that
  fallback manufactured **four false microvariants in 202 called alleles**
  (TH01 6.3, TH01 6.1, D13S317 14.1, D18S51 11.3 — all integers under real POA).
  The fallback still exists so the package stays importable, but it raises a
  `CONSENSUS_FALLBACK` warning on every affected marker.
- **Reads are binned by repeat-core length, not by window length.** Panel
  windows are ~80% flank, so binning on raw length let flank indel errors split
  one allele into two clusters. Core binning cut cluster fragmentation by 29%
  and eliminated all false mixture flags on the test set. Binning by repeat
  *unit count* was evaluated and rejected — `[AATG]9` and `[AATG]6 ATG [AATG]3`
  are both nine units, so it would have destroyed TH01 9 vs 9.3.
- **Stutter is modelled log-linearly in LUS**, fitted to ONT data
  (R² 0.965), not as a flat per-marker rate. The measured rate spans more than
  10× across LUS 10–14, which no constant can express. See
  [`docs/stutter_calibration.md`](docs/stutter_calibration.md) — including the
  caveat that the shipped model is fitted on **PCR-free WGS** and must be
  re-fitted before amplicon casework.
- **Haplotype tags are used as evidence.** On a phased BAM, two candidates
  confidently on the same haplotype and within 2 bp of each other cannot both
  be real alleles, so the weaker is suppressed and recorded.
- **Every QC condition that is checked is reported as checked.** The audit
  record lists every flag code the pipeline can raise, so a code absent from
  the counts means "evaluated and not found" rather than "never looked at".

## Requirements

- Python ≥ 3.11
- A **POA backend** — `pyspoa` by default. `pyabpoa` is also supported but does
  not build on macOS arm64 (it hardcodes AVX2 x86 intrinsics).
- `samtools` for preparing inputs; `LongTR` only if you want Layer 1.

```bash
pip install -e '.[poa]'
```

Then confirm the backend resolved, because a run without one is a different run:

```bash
frontstr doctor --bam sample.bam --panel examples/panels/codis_20_grch38.yaml
```

## Usage

FRONTStr currently takes an **indexed, aligned BAM or CRAM**. Align your reads
first (`minimap2 -ax map-ont`); see [Limitations](#limitations).

Call a profile and print it:

```bash
frontstr interpret --bam sample.bam --panel examples/panels/codis_20_grch38.yaml
```

Full pipeline with exports, an HTML report and an audit trail:

```bash
frontstr export --bam sample.bam --panel examples/panels/codis_20_grch38.yaml --out-dir out/ --formats profile,evidence,seqs,json,html --operator "your name" --run-id RUN-001
```

That writes `<sample>.profile.csv`, `<sample>.evidence.csv`,
`<sample>.seqs.csv`, `<sample>.json`, a self-contained `<sample>.html`, and a
`frontstr.log.jsonl` process log. The canonical JSON carries a sealed audit
record: tool version, POA backend, stutter model, every threshold that moved a
call, input hashes and the flag census.

### Exports

`--formats` accepts `profile`, `evidence`, `seqs` (CSV), `json`,
`json-compact`, `html`, `xlsx` and `vcf`.

- **`vcf`** — a native, sequence-resolved VCF. `ALT` is the allele's sequence,
  not a length, so iso-alleles remain distinct records; the repeat count rides
  along in `FORMAT/MC` and FRONTStr's integer per-allele coverage in
  `FORMAT/AD`. QC warnings become `FILTER` entries. It is bgzip/tabix
  indexable and queryable with bcftools, which is the point — it exists so
  FRONTStr can be benchmarked against other callers. **Needs `--reference`.**
- **`xlsx`** — a five-sheet review workbook (Profile, Sequences, Evidence, QC,
  Audit). Markers carrying a warning are tinted, QC is ordered worst-first, and
  the Evidence sheet keeps the clusters that were *not* called.

Other commands: `evidence` (per-locus cluster dump for debugging), `call`
(LongTR), `batch` (multi-sample from a manifest), `calibrate-stutter`
(fit a stutter model to your own data), `doctor`, `inspect`.

Run `frontstr --help` for the full list.

## Limitations

Read this section as part of the documentation, not as a disclaimer.

- **`frontstr run` is a stub.** The FASTQ → alignment → report path is not
  wired; `ingest.align` raises. Align externally and start from a BAM.
- **Not laboratory-validated.** No forensic partner has signed off on a
  validation report. There is no mixture series, no dropout study and no NIST
  control run. The reference profile it is regression-tested against comes from
  HipSTR on matched Illumina data — that is caller-vs-caller concordance, not
  validation against a reference method.
- **The stutter model is PCR-free.** Fitted on WGS data, so it contains no PCR
  slippage component and will under-predict stutter on an amplicon panel. The
  model records its own protocol so this cannot be lost.
- **vWA and D21S11 do not report the legacy kit allele number.** FRONTStr
  reports the sequence-derived repeat count; for these two compound loci that
  is not the CE-kit designation, and no single correction reconciles them. Both
  raise a `CE_NOMENCLATURE_OFFSET` warning with the reason. Do not compare
  those numbers directly against a CE profile.
- **The analytical and calling thresholds (0.02 / 0.10) are not data-derived.**
  They are chosen defaults. The low-coverage floor and the stutter model *are*
  derived; these two are not, and that is a gap.
- **No CODIS `.cmf`, MIDST or PDF export**, and no tidy/Parquet dataset for
  cohort-scale analysis, despite what earlier plans promised.
- **Illumina data will not work** against the shipped panel. Its windows are
  ±100 bp and the pileup requires reads that fully span them; 150 bp reads
  cannot. This is a deliberate ONT-first design choice, not an oversight.

See [`ROADMAP.md`](ROADMAP.md) for what is planned and what is deliberately
deferred.

## Development

```bash
pip install -e '.[dev,poa]'
pytest                      # unit suite
pytest -m integration       # end-to-end regression, needs local ONT slices
ruff check frontstr tests && ruff format --check frontstr tests && mypy frontstr
```

The integration tests run the whole pipeline against a real ONT R10 BAM and
assert the genotype of every marker. The slices are not versioned (~200 MB), so
those tests skip cleanly when the data is absent.

## License

MIT — see [`LICENSE`](LICENSE).

## Citing

If you use FRONTStr in research, please cite:

> [TBD — pending preprint]

### Prior art

FRONTStr is an independent implementation. It uses and draws on the following:

- **LongTR** — Tang et al. Used directly as the optional Layer 1 genotyper.
- **ISFG sequence nomenclature** — the bracketed-repeat convention FRONTStr
  emits follows the ISFG recommendations and the counting rules described by
  Phillips (2018).
- **toaSTR** — Ganschow S, Silvery J, Kalinowski J, Tiemann C (2018), *toaSTR:
  A web application for forensic STR genotyping by massively parallel
  sequencing*, Forensic Sci Int Genet 37:21–28. Cited as a methodological
  reference for sequence-based STR genotyping, in particular the LUS/SLUS
  framing of stutter. **FRONTStr is not a successor to, port of, or derivative
  of toaSTR, and contains none of its code** — toaSTR's public repository
  distributes only a Docker Compose file and a SQL schema, and no toaSTR source
  is present in this project. The stutter rates FRONTStr once carried over from
  the CE/short-read literature have since been replaced with values fitted to
  ONT data.
- **STRspy** — Hall, Kesharwani et al. Referenced for allele-database design.
