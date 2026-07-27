# FRONTStr

> **Forensic Ranked Output for Nanopore Tandem Short Tandem Repeats**

A forensic STR profiling toolkit for **long-read sequencing**, built ONT-first:
every default, calibration and threshold in it was measured on Oxford Nanopore
R10 data rather than carried over from capillary electrophoresis or short-read
practice.

That is the point of the project. A long read spans an entire STR locus and its
flanks at once, which makes it possible to report the **sequence** of each
allele rather than only its length — iso-alleles, microvariants and interrupted
repeat structures that CE cannot resolve. But a nanopore read's error profile
is nothing like a capillary trace's, and a caller that inherits CE-era
constants is confidently wrong in ways that are hard to see. Where FRONTStr had
inherited such constants, they have been replaced with values fitted to real
ONT data, and the measurements are in [`docs/`](docs/) — including the ones
that changed the answer.

## Status

**Pre-alpha.** The pipeline works end to end from an aligned BAM and is
regression-tested against a real ONT sample. It has **not** been validated by a
forensic laboratory and is not fit for casework. The
[Limitations](#limitations) are specific and worth reading before you draw
conclusions from any output.

## What it does

**FRONTStr calls genotypes by itself.** Every allele, read count and sequence
in its output comes from its own pileup of the BAM.

| Stage | Module | What it does |
|---|---|---|
| Evidence | `frontstr.evidence` | Per-locus pileup → clustering by repeat-core length → POA consensus. **Integer per-allele read counts straight from the BAM.** |
| Interpretation | `frontstr.interp` | ISFG nomenclature, allele numbering, stutter, haplotype-aware suppression, genotype calling. |
| QC & audit | `frontstr.interp.qc`, `frontstr.audit` | Coverage, strand-bias and nomenclature flags; the sealed audit record. |

Coverage is counted, not divided out of a caller's `BPDIFFS` field, and the
ISFG bracket string is computed from the cluster consensus rather than from a
VCF `ALT`.

> **LongTR is optional and off by default.** `frontstr.caller` wraps it, but it
> runs only if you pass `--longtr-vcf`, and all it does then is cross-check —
> disagreement raises a flag for an analyst, never changes a call. FRONTStr
> does not need LongTR installed.

Four choices drive most of what the tool reports, each measured rather than
asserted:

- **POA consensus is required, not optional.** Without it the consensus is a
  single unpolished read, which on the test set manufactured four false
  microvariants in 202 called alleles.
- **Reads are binned by repeat-core length**, not window length — panel windows
  are ~80% flank. Binning by repeat *unit count* was evaluated and rejected: it
  would merge TH01 9 with 9.3.
- **Stutter is log-linear in LUS** (R² 0.965), because the measured rate spans
  more than 10× across LUS 10–14.
- **Haplotype tags are evidence.** Two candidates on the same haplotype within
  2 bp cannot both be real alleles.

**→ [`docs/how_it_works.md`](docs/how_it_works.md) walks every stage with
worked examples, the derivation behind each threshold, and the full current
state of the project.**

## Install

Requires Python ≥ 3.11 and a POA backend.

```bash
pip install -e '.[poa]'
```

`pyspoa` is the default backend; `pyabpoa` also works but does not build on
macOS arm64 (it hardcodes AVX2 x86 intrinsics). Confirm the backend resolved
before trusting a run — one without it is a different run:

```bash
frontstr doctor --bam sample.bam --panel examples/panels/codis_20_grch38.yaml
```

## Usage

FRONTStr takes an **indexed, aligned BAM or CRAM**. Align first with
`minimap2 -ax map-ont` — the FASTQ path is not wired yet.

Call a profile and print it:

```bash
frontstr interpret --bam sample.bam --panel examples/panels/codis_20_grch38.yaml
```

Full pipeline with exports and an audit trail:

```bash
frontstr export --bam sample.bam --panel examples/panels/codis_20_grch38.yaml --out-dir out/ --formats profile,seqs,json,html,xlsx,vcf --reference GRCh38.fa --operator "your name" --run-id RUN-001
```

Alongside the exports this writes a `frontstr.log.jsonl` process log, and the
canonical JSON carries a sealed audit record: tool version, POA backend,
stutter model, every threshold that moved a call, input hashes and the flag
census.

### Export formats

| `--formats` | Output |
|---|---|
| `profile` | CSV, one row per marker — the genotype table |
| `evidence` | CSV, one row per cluster, including the uncalled ones |
| `seqs` | CSV, one row per called allele with ISFG and consensus |
| `json` / `json-compact` | The canonical record: everything, plus the audit block |
| `html` | Self-contained offline report |
| `xlsx` | Five-sheet review workbook (Profile, Sequences, Evidence, QC, Audit) |
| `vcf` | Native sequence-resolved VCF — **needs `--reference`** |

The VCF is native, not an annotation of LongTR's output. `ALT` carries the
allele's sequence so iso-alleles stay distinct; the repeat count rides in
`FORMAT/MC` and per-allele coverage in `FORMAT/AD`. It is bgzip/tabix
indexable and bcftools-queryable, which is the point — it exists so FRONTStr
can be benchmarked against other callers.

Other commands: `evidence` (per-locus cluster dump), `batch` (multi-sample from
a manifest), `calibrate-stutter` (fit a stutter model to your own data), `call`
(LongTR), `doctor`, `inspect`. Run `frontstr --help` for the full list.

## Limitations

Part of the documentation, not a disclaimer.

- **`frontstr run` is a stub.** The FASTQ → alignment → report path is not
  wired; `ingest.align` raises. Align externally and start from a BAM.
- **Not laboratory-validated.** No forensic partner has signed off. No mixture
  series, no dropout study, no NIST control. The reference profile it is tested
  against comes from HipSTR on matched Illumina data — caller-vs-caller
  concordance, not validation against a reference method.
- **The stutter model is PCR-free.** Fitted on WGS, so it has no PCR slippage
  component and will under-predict stutter on an amplicon panel. The model
  records its own protocol so this cannot be lost.
- **vWA and D21S11 do not report the legacy kit allele number.** FRONTStr
  reports the sequence-derived repeat count; for these two compound loci that
  is not the CE-kit designation, and no single correction reconciles them. Both
  raise `CE_NOMENCLATURE_OFFSET`. Do not compare those numbers against a CE
  profile directly.
- **The analytical and calling thresholds (0.02 / 0.10) are not data-derived.**
  The coverage floor and the stutter model are; these two are chosen defaults,
  and that is a gap.
- **No CODIS `.cmf`, MIDST or PDF export**, and no tidy/Parquet dataset for
  cohort-scale analysis.
- **Illumina data will not work** against the shipped panel: its windows are
  ±100 bp and the pileup needs reads that span them fully. A deliberate
  ONT-first choice, not an oversight.

`docs/how_it_works.md` §12 lists the known gaps *inside* what does work.
[`ROADMAP.md`](ROADMAP.md) covers what is planned and what is deferred.

## Documentation

| Document | Contents |
|---|---|
| [`docs/how_it_works.md`](docs/how_it_works.md) | Every stage step by step, with examples and measurements; full current state |
| [`docs/stutter_calibration.md`](docs/stutter_calibration.md) | How the stutter model was measured, and its caveats |
| [`docs/architecture.md`](docs/architecture.md) | Module layout |
| [`ROADMAP.md`](ROADMAP.md) | Delivery plan |

## Development

```bash
pip install -e '.[dev,poa]'
pytest                      # unit suite
pytest -m integration       # end-to-end regression, needs local ONT slices
ruff check frontstr tests && ruff format --check frontstr tests && mypy frontstr
```

The integration tests run the whole pipeline against a real ONT R10 BAM and
assert the genotype of every marker. Those slices are not versioned (~200 MB),
so the tests skip cleanly when the data is absent.

## License

MIT — see [`LICENSE`](LICENSE).

## Citing

If you use FRONTStr in research, please cite:

> [TBD — pending preprint]

### Prior art

FRONTStr is an independent implementation. It uses and draws on:

- **LongTR** — Tang et al. Wrapped as an optional cross-check.
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
  is present in this project. The stutter rates FRONTStr once carried from the
  CE/short-read literature have been replaced with values fitted to ONT data.
- **STRspy** — Hall, Kesharwani et al. Referenced for allele-database design.
