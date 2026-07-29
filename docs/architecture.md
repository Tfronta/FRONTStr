# Architecture

The module map and the dependency rules between them. For what the pipeline
*does* at each stage, see [`how_it_works.md`](how_it_works.md); this document
is about where the code lives and why.

## The dependency rule

Dependencies run **one way**: Evidence may not import Interpretation, and
neither may import the report or export layers. The rule is enforced by review,
not by tooling, but it is currently unbroken — `frontstr.evidence` and
`frontstr.caller` contain no import of `frontstr.interp`.

It buys two things. Changing a forensic threshold touches only Interpretation.
Changing how reads are extracted or clustered touches only Evidence, and cannot
quietly alter a calling rule.

Two modules exist specifically to keep that rule intact:

- **`frontstr/motifs.py`** — motif-run scanning, repeat-core location,
  reverse-complement. Pure sequence functions with no forensic opinion, sitting
  outside the layering because *all three* of Evidence, Interpretation and
  Panel need them. Evidence bins reads by repeat-core length; Interpretation
  needs the same primitives for stutter, ISFG and catalog lookup. Without this
  module, Evidence would have to import Interpretation.
- **`frontstr/log.py`** — structured logging, deliberately free of any FRONTStr
  import. Logging and the audit record started as one module and produced a
  circular import: `interp/profile.py` needs a logger, while the audit record
  needs the interp models. Splitting them is both the fix and the correct
  layering — the domain should not depend on the audit trail.

---

## Evidence — `frontstr.evidence`

Reads to allele candidates. Everything here is measurement; nothing decides
whether something is an allele.

| Module | Responsibility |
|---|---|
| `pileup` | Per-locus extraction with pysam. One `Observation` per read that spans the whole window, carrying sequence in reference orientation, `HP` tag, strand and mean quality. |
| `cluster` | Two-stage grouping: bin by **repeat-core length**, then merge by edit-distance identity within the bin. |
| `consensus` | POA backend chain (`pyabpoa` → `pyspoa` → mode fallback) using **global** alignment. Records which method produced each consensus. |

Output: `Cluster` objects with integer `n_reads`, per-haplotype counts, a
consensus sequence and a `consensus_method`.

Global alignment is not an incidental choice. Clusters are length-binned, so
members are equal-length full-window sequences; local alignment would be free
to trim flanks and silently change the called length, which *is* the CE number.

---

## Interpretation — `frontstr.interp`

Allele candidates to a forensic genotype. Every decision rule lives here.

| Module | Responsibility |
|---|---|
| `models` | `Allele`, `MarkerResult`, `Flag`, `IsoAllele` and the enums. The canonical allele number, its derivation method and its display label are computed here so no view can disagree about them. |
| `isfg` | Bracketed nomenclature compression; the legacy CE from length or bracket count, now a fallback. |
| `naming` | STRNaming-backed allele naming — the canonical CE. Offline, from a committed GRCh38 slice cache. |
| `allele_numeric` | Reference-anchored allele numbering for compound markers. |
| `stutter` | Builds expected stutter coverage per virtual stutter sequence, from a calibrated `StutterModel`. |
| `classify` | `allele \| stutter \| artefact \| noise \| inexact_allele \| deletion \| hp_phantom`. |
| `haplotype` | Suppresses same-haplotype split-allele phantoms. No-op on unphased BAMs. |
| `catalog` | Annotates alleles against a curated iso-allele catalog (optional). |
| `triallelic` | `call_profile` — 1, 2 or 3+ alleles, plus mixture suspicion. |
| `amel` | Amelogenin sex typing. Not a tandem repeat, so it bypasses the STR path entirely. |
| `flags` | Marker-level flags intrinsic to a finished call (triallelic, iso-allele, unpolished consensus). |
| `qc` | Run-level flags that depend on a laboratory threshold (coverage, strand bias, kit nomenclature). |
| `profile` | The orchestrator: `interpret_marker` and `interpret_run`. |

---

## Caller — `frontstr.caller` (present, not wired)

Wraps LongTR. **Nothing in the pipeline calls it and no command exposes it.**
It was once an optional cross-check behind `--longtr-vcf`; that wiring, and the
`interp.concordance` glue it fed, were removed. Benchmarking FRONTStr against
another caller is a separate exercise from running FRONTStr, and doing both at
once made it hard to say which tool produced a number.

The code stays because the argv construction, BED emission and VCF parsing are
the tedious parts and are worth keeping for a future benchmark harness. It is
importable and tested; it is simply not cabled to anything.

| Module | Responsibility |
|---|---|
| `longtr` | argv construction and subprocess invocation |
| `vcf` | cyvcf2 → in-memory result objects |

What this layer never produces: per-allele coverage. That comes from Evidence.

---

## Cross-cutting

| Module | Responsibility |
|---|---|
| `panel.models`, `panel.loader` | Versioned panel and marker definitions from YAML |
| `panel.catalog` | Iso-allele catalog model and JSON I/O |
| `panel.calibrate` | Derives per-marker `corr_value` from a reference FASTA |
| `panel.stutter_calib` | Measures stutter from real BAMs and fits a `StutterModel` |
| `panel.seed_strseq` | STRSeq catalog builder — assembly done, the NCBI fetch is a stub |
| `panel.bed` | Panel → BED. Lives here, not in `caller`: this is how a panel serializes, and the report embeds it |
| `ingest.detect`, `ingest.validate` | Input format sniffing and BAM header checks |
| `ingest.align` | minimap2 wrapper — **not implemented**, raises |
| `report.payload` | `serialize_run`: the single serializer every consumer shares |
| `report.html`, `svg_charts`, `ngs_display` | Self-contained HTML report |
| `exports.csv`, `json`, `xlsx`, `vcf` | Output formats |
| `audit` | `AuditRecord`: run configuration, flag census, integrity seal |
| `params` | Every run knob, its default and that default's provenance (`derived` / `chosen` / `convention`). Top level: the knobs span Evidence and Interpretation, and both the CLI echo and the report's parameter table read from it |
| `log` | JSONL process log; renders JSONL to file and readable prose to a terminal |
| `trace` | Per-locus narrative behind `interpret --trace`. `LocusTrace` is a plain record the evidence and interp layers fill; `render_locus` turns it into prose. Top level because it spans both layers, and split record-from-rendering so the same trace can later feed the audit record or the HTML report |
| `batch` | Multi-sample orchestration from a manifest |
| `cli` | Typer entry point; every subcommand |
| `errors` | The `FrontstrError` hierarchy |
| `motifs` | Motif-run primitives shared by all layers — see the dependency rule above |

Not built, despite appearing in earlier plans: PDF report, batch report
generator, CODIS CMF, NIST MIDST, ZIP bundle, tidy/Parquet dataset.

---

## Data flow

```
FASTQ  ──►  ingest.align_to_reference        ✗ NOT IMPLEMENTED — align externally
                                                (minimap2 -ax map-ont)
indexed BAM / CRAM
        │   single source of truth for everything downstream
        │
        └──► evidence.pileup
                 │
                 ▼
             evidence.cluster
                 │
                 ▼
             evidence.consensus
                 │
                 ▼
             interp.profile
               naming → isfg → stutter → classify
               → haplotype → catalog → triallelic
                 │
                 ▼
             interp.flags  →  interp.qc
                 │
                 ▼
             MarkerResult (per locus) + Allele (per cluster)
                 │
                 ▼
             report.payload.serialize_run  +  audit.build_audit_record
                 │
                 ▼
             report.html  ·  exports.{csv, json, xlsx, vcf}  ·  log.jsonl
```

---

## Why the layering matters forensically

- **Provenance is traceable.** Every cluster knows the read IDs that produced
  it, every allele knows its cluster, and every run records the software
  versions, backends and thresholds that produced it in a sealed audit record.
- **Independent of any external caller.** Coverage comes from Evidence, and no
  external caller is wired in at all, so there is nothing whose version or
  absence could move the numbers.
- **Thresholds are defensible and located.** Stutter rates live in a versioned
  `StutterModel`, QC thresholds in a `QcThresholds` object, marker parameters
  in the panel YAML. None of them are buried in caller flags, and all of them
  are serialized into the audit record of every run that used them.
