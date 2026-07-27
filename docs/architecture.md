# Architecture

FRONTStr is organised as **three layers with strict, one-way dependencies**.
The boundaries are not negotiable: if LongTR ships a new version, only Layer 1
needs work; if the lab changes its forensic thresholds, only Layer 3 changes.

## Layer 1 — Caller (`frontstr.caller`)

Wraps the external caller (LongTR) and parses its VCF.

- `caller.longtr.build_longtr_argv` — deterministic argv for a given run.
- `caller.longtr.run_longtr` — subprocess invocation.
- `caller.vcf` — cyvcf2 → in-memory result objects.

What this layer **never** produces: per-allele integer coverage. That is the
job of Layer 2.

## Layer 2 — Evidence (`frontstr.evidence`)

Sequence-level pileup straight from the BAM.

- `evidence.pileup.pileup_locus` — pysam-driven per-locus extractor.
- `evidence.cluster.cluster_observations` — length-binned + edit-distance
  clustering with POA local consensus.

Result: a list of `Cluster` objects per locus with integer `n_reads`,
per-haplotype counts (`HP1`/`HP2`/none), and a consensus sequence ready for
ISFG compression.

## Layer 3 — Interpretation (`frontstr.interp`)

Forensic decision rules.

- `interp.isfg` — bracketed nomenclature compression.
- `interp.stutter` — LUS/SLUS expected-stutter model + per-marker overrides.
- `interp.classify` — `allele | stutter | artefact | noise | inexact | deletion`.
- `interp.triallelic` — `call_profile` supporting 1, 2, or 3+ alleles plus
  mixture detection.
- `interp.concordance` — cross-check with LongTR's GT; flag discordances.

## Cross-cutting

- `panel` — versioned forensic panel + (later) allele catalog.
- `ingest` — input format detection, BAM validation, alignment.
- `report` — HTML / PDF / batch generators sharing one serializer.
- `exports` — CSV, XLSX, VCF-extended, CODIS CMF, MIDST, ZIP bundle.

## Data flow

```
input file (fastq | bam | cram)
        │
        ▼
ingest.detect_input → ingest.validate_bam OR ingest.align_to_reference
        │
        ▼
sorted, indexed BAM (single source of truth for everything downstream)
        │
        ├──► caller.longtr  ─────────────►  tr_calls.vcf
        │
        └──► evidence.pileup
                 │
                 ▼
             cluster_observations
                 │
                 ▼
             interp.isfg + interp.classify + interp.triallelic
                 │
                 ▼
             Result(per locus) + Allele(per cluster)
                 │
                 ▼
             report.html  +  exports.{csv,xlsx,vcf_extend,codis,bundle}
```

## Why this layering matters forensically

- **Auditable provenance**: every cluster knows the read IDs that produced it,
  every allele knows which cluster it came from, every result knows which
  pipeline parameters and software versions produced it.
- **Independent of caller**: if LongTR shifts, our coverage numbers don't
  change. Layer 2 is the source of truth.
- **Defensible thresholds**: stutter rules live in Layer 3 and are versioned
  with the panel, not buried in caller flags.
