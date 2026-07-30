# Architecture

Where the code lives. For what the pipeline does at each stage, see
[`how_it_works.md`](how_it_works.md).

## The dependency rule

Dependencies run one way. Evidence does not import Interpretation, and neither
imports the report or export layers. The rule is enforced by review rather than
tooling.

Changing a forensic threshold touches only Interpretation. Changing how reads
are extracted or clustered touches only Evidence, and cannot alter a calling
rule.

Two modules sit outside the layering to keep it intact:

- **`frontstr/motifs.py`**: motif-run scanning, repeat-core location,
  reverse-complement. Pure sequence functions with no forensic decision in
  them. Evidence, Interpretation and Panel all need these primitives; without a
  shared module Evidence would have to import Interpretation.
- **`frontstr/log.py`**: structured logging, with no FRONTStr import of its
  own. The domain does not depend on the audit trail.

---

## Evidence: `frontstr.evidence`

Reads to allele candidates. Measurement only; nothing here decides whether
something is an allele.

| Module | Responsibility |
|---|---|
| `pileup` | Per-locus extraction with pysam. One `Observation` per read spanning the whole window, carrying sequence in reference orientation, `HP`/`PS` tags, strand and mean quality |
| `cluster` | Two-stage grouping: bin by repeat-core length, then split by edit-distance identity within the bin |
| `consensus` | POA backend chain (`pyspoa`, `pyabpoa`, mode fallback), global alignment. Records which method produced each consensus |

Output: `Cluster` objects with integer `n_reads`, per-haplotype counts, a
consensus sequence and a `consensus_method`.

Alignment is global. Clusters are length-binned, so members are equal-length
full-window sequences, and local alignment could trim flanks and change the
called length, which is the allele number.

---

## Interpretation: `frontstr.interp`

Allele candidates to a forensic genotype. Every decision rule lives here.

| Module | Responsibility |
|---|---|
| `models` | `Allele`, `MarkerResult`, `Flag`, `IsoAllele` and the enums. The canonical allele number, its derivation method and its display label are computed here, so no view can disagree about them |
| `naming` | STRNaming-backed allele naming, offline from a committed GRCh38 slice cache |
| `isfg` | Bracketed nomenclature compression, and the fallback allele number from length or bracket count |
| `allele_numeric` | Reference-anchored allele numbering for compound markers |
| `stutter` | Expected stutter coverage per virtual stutter sequence, from a calibrated `StutterModel` |
| `classify` | `allele`, `stutter`, `artefact`, `noise`, `inexact_allele`, `deletion`, `hp_phantom` |
| `haplotype` | Same-haplotype phantom suppression, and the opposite-haplotype test. No-op on unphased BAMs |
| `catalog` | Annotates alleles against a curated iso-allele catalog (optional) |
| `triallelic` | `call_profile`: one, two, or three or more alleles, plus mixture suspicion |
| `amel` | Amelogenin sex typing. Not a tandem repeat, so it bypasses the STR path |
| `flags` | Marker-level flags intrinsic to a finished call |
| `qc` | Run-level flags that depend on a laboratory threshold: coverage, strand bias, allele balance |
| `profile` | Orchestrator: `interpret_marker` and `interpret_run` |

---

## Cross-cutting

| Module | Responsibility |
|---|---|
| `panel.models`, `panel.loader` | Versioned panel and marker definitions from YAML |
| `panel.bed` | Panel to BED, and BED to panel for `--bed` |
| `panel.catalog` | Iso-allele catalog model and JSON I/O |
| `panel.calibrate` | Derives per-marker `corr_value` from a reference FASTA |
| `panel.stutter_calib` | Measures stutter from BAMs and fits a `StutterModel` |
| `panel.seed_strnaming` | Builds the committed GRCh38 slice cache used by `interp.naming` |
| `panel.seed_strseq` | STRSeq catalog builder. Assembly done, the NCBI fetch is a stub |
| `ingest.detect`, `ingest.validate` | Input format sniffing and BAM header checks |
| `ingest.align` | minimap2 wrapper. Not implemented, raises |
| `report.payload` | `serialize_run`, the single serializer every consumer shares |
| `report.html`, `svg_charts`, `ngs_display` | Self-contained single-sample HTML report |
| `report.cohort` | Multi-sample view: one block per marker, one row per sample |
| `exports.csv`, `json`, `xlsx`, `vcf` | Output formats |
| `exports.tidy` | Cohort-scale long table, CSV and Parquet, built from run JSONs |
| `audit` | `AuditRecord`: run configuration, flag census, integrity seal |
| `params` | Every run knob, its default and that default's provenance. Top level because the knobs span Evidence and Interpretation |
| `log` | JSONL process log, rendered as JSONL to file and readable text to a terminal |
| `trace` | Per-locus narrative behind `--trace`. `LocusTrace` is a plain record the evidence and interp layers fill; `render_locus` turns it into text. Top level because it spans both layers |
| `batch` | Multi-sample orchestration from a manifest |
| `cli` | Typer entry point and every subcommand |
| `errors` | The `FrontstrError` hierarchy |
| `motifs` | Motif-run primitives shared by all layers |

`frontstr.caller` wraps an external caller. Nothing in the pipeline calls it and
no command exposes it.

Not built: PDF report, CODIS CMF, NIST MIDST.

---

## Data flow

```
FASTQ  ──►  ingest.align_to_reference        NOT IMPLEMENTED, align externally
                                             (minimap2 -ax map-ont)
indexed BAM / CRAM
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
             report.html · report.cohort · exports.{csv, json, xlsx, vcf, tidy}
```

---

## Provenance

Every cluster records the read IDs that produced it, every allele records its
cluster, and every run records the software versions, backends and thresholds
that produced it in the audit record.

Stutter rates live in a versioned `StutterModel`, QC thresholds in a
`QcThresholds` object, and marker parameters in the panel YAML. All of them are
serialized into the audit record of the run that used them.
