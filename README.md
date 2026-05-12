# FRONTStr

> **Forensic Ranked Output for Nanopore Tandem Short Tandem Repeats**

FRONTStr is a forensic STR profiling toolkit specifically designed for long‑read
sequencing data (Oxford Nanopore R10, PacBio HiFi). It is the successor to
[toaSTR](https://www.toastr.app/) for the long‑read era, with a stricter
forensic mindset:

- **True per‑allele coverage** — integer read counts from a sequence‑level
  pileup, not approximations from `BPDIFFS`/`PDP` divisions.
- **ISFG bracketed nomenclature** on the cluster consensus, not on the VCF ALT.
- **Tri‑allelic patterns supported out of the box** (TPOX type I/II, vWA, FGA, …).
- **Self‑contained, visually rich HTML reports** — single‑file, offline, emailable.
- **Chain‑of‑custody by design** — every artefact hashed; deterministic POA seeds.
- **CODIS `.cmf` export, NIST MIDST**, multi‑sheet XLSX, full case ZIP bundle.

## Status

Pre‑alpha. See [`ROADMAP.md`](ROADMAP.md) for phased delivery plan.

## Quick start (once installed)

```bash
# Single sample (ONT FASTQ → full pipeline → HTML report)
frontstr run \
    --sample S001 \
    --platform ont \
    --panel examples/panels/codis_20_grch38.yaml \
    --reference /refs/GRCh38.p14.fa \
    --input sample.fastq.gz \
    --out  out/S001/

# Single sample (pre‑aligned BAM, passthrough)
frontstr run \
    --sample S001 --platform ont \
    --panel examples/panels/codis_20_grch38.yaml \
    --reference /refs/GRCh38.p14.fa \
    --input sample.sorted.bam \
    --out out/S001/

# Batch (12+ samples from one MinION flowcell)
frontstr batch \
    --manifest batch.tsv \
    --panel examples/panels/codis_20_grch38.yaml \
    --reference /refs/GRCh38.p14.fa \
    --out out/batch-2026-05-11/
```

## Architecture

Three layers with strict responsibilities (see [`docs/architecture.md`](docs/architecture.md)):

```
┌─────────────────────────────────────────────────────────┐
│  Layer 1 — Caller (LongTR): VCF for QC                  │
├─────────────────────────────────────────────────────────┤
│  Layer 2 — Evidence (sequence pileup): real coverage    │
├─────────────────────────────────────────────────────────┤
│  Layer 3 — Interpretation (ISFG, tri‑allelic, mixture)  │
└─────────────────────────────────────────────────────────┘
```

## Why FRONTStr exists

| Need | toaSTR v1 | LongTR | STRspy 2.0 | FRONTStr |
|---|:--:|:--:|:--:|:--:|
| Long‑read native | ❌ | ✅ | ✅ | ✅ |
| ISFG nomenclature | ✅ Illumina | ❌ | ❌ | ✅ |
| Per‑allele integer coverage | ✅ Illumina | ⚠ approx | ⚠ raw counts | ✅ |
| Tri‑allelic patterns | ❌ | ❌ | ❌ (top‑2) | ✅ |
| Self‑contained HTML report | ❌ | ❌ | ❌ | ✅ |
| CODIS `.cmf` export | ❌ | ❌ | ❌ | ✅ |
| Chain‑of‑custody hashes | ❌ | ❌ | ❌ | ✅ |
| Mixture detection (batch) | ❌ | ❌ | ❌ | ✅ |

## License

MIT (TBD; pending review with stakeholders).

## Citing

If you use FRONTStr in casework or research, please cite:

> [TBD — pending preprint]

Methodologically FRONTStr builds on:

- Tang et al., *LongTR* (2024) — long‑read TR genotyper, used as Layer 1.
- Ganschow et al., *toaSTR* (2018) — forensic Illumina STR; LUS/SLUS stutter model.
- Hall, Kesharwani et al., *STRspy 2.0* (2026) — allele‑DB design inspiration.
