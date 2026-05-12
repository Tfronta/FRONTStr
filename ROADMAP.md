# FRONTStr — Roadmap

> Roadmap from "we have a plan" to "v1.0 production".
> See [`../plan-longtr-improved.md`](../plan-longtr-improved.md) for the full design.

## Strategy: ship the killer feature first

The unique forensic value of FRONTStr is **per‑allele integer coverage + ISFG +
beautiful HTML report**. Everything else (web UI, batch, mixture detection, CODIS
export, etc.) is high value but **strictly secondary**. We deliver the killer
feature **first**, get one forensic lab to validate it, then add the rest.

Concretely: **no FastAPI, no React, no Postgres, no Redis, no Docker
orchestration in the MVP**. CLI + HTML output. That is enough to validate the
science with real labs in week 4.

---

## Phase 0 — Bootstrap (week 1, 3–5 days)

**Goal**: empty repo → running test suite.

- [x] Directory structure: `frontstr/`, `tests/`, `docker/`, `examples/`, `docs/`
- [x] `pyproject.toml` with pinned deps (Python 3.11+)
- [x] `.gitignore`, `LICENSE`, `README.md`
- [ ] `conftest.py` with shared fixtures (tmp paths, panels, etc.)
- [ ] `pre‑commit` config: `ruff`, `black`, `mypy --strict`
- [ ] GitHub Actions: lint + test on push
- [ ] `frontstr/__main__.py` runnable: `python -m frontstr --help`

**Acceptance**: `pytest` passes (with no tests yet); `frontstr --help` prints.

---

## Phase 1 — MVP CLI single‑sample (weeks 1–3)

**Goal**: `frontstr run -i sample.fastq.gz -p panel.yaml --reference ref.fa
--out out/` produces a working `report.html` with the unique value props.

### 1.1 Ingest (3 days)

- [ ] `ingest/detect.py` — magic‑bytes sniff (FASTQ / BAM / CRAM / uBAM)
- [ ] `ingest/validate.py` — BAM header sanity (@RG, @SQ MD5, MAPQ distribution)
- [ ] `ingest/align.py` — minimap2 wrapper (`-ax map-ont` / `map-hifi`)
- [ ] Tests with synthetic BAM/FASTQ fixtures
- [ ] CLI subcommand `frontstr ingest` (dev convenience)

### 1.2 Panel model (1 day)

- [ ] `panel/models.py` — `Panel`, `System`, `AlleleCatalogEntry` dataclasses
- [ ] `panel/loader.py` — YAML → in‑memory model
- [ ] `examples/panels/codis_20_grch38.yaml` — 24 markers seeded
- [ ] Tests: load, validate, dump

### 1.3 Caller (LongTR) — DONE

- [x] `caller/longtr.py` — argv builder + subprocess runner + `LongTRRunner` orchestrator
- [x] `caller/bed.py` — single panel BED and per-chromosome split
- [x] `caller/vcf.py` — cyvcf2 parsing into `LongTRResult` / `LongTRAlleleSpec` / `LongTRSampleCall`
- [x] CLI `frontstr call` with `--parse-only` for offline VCF inspection
- [x] Unit tests for argv (ONT vs HiFi), BED writer, VCF parser (het, hom, `<DEL>`, INEXACT_ALLELE, phased, missing GT), subprocess failure paths

### 1.4 Evidence layer (the unique value) (4–5 days)

- [ ] `evidence/pileup.py` — `pileup_locus` with pysam
- [ ] `evidence/cluster.py` — cluster by length + edit distance + POA consensus
  (use `pyabpoa`; fallback to `spoa` if not installable)
- [ ] Tests with synthetic ONT reads: known clusters in/out
- [ ] CLI subcommand `frontstr evidence <bam> <bed>` (debug)

### 1.5 Interpretation — DONE

- [x] `interp/models.py` — `Allele`, `MarkerResult`, `AlleleStatus`, `CallRule`, `TriType`
- [x] `interp/isfg.py` — `compress_isfg` + `ce_from_length`
- [x] `interp/stutter.py` — LUS/SLUS detection + `build_expected_stutter` (with overrides)
- [x] `interp/classify.py` — coverage + sequence-based two-rule classifier
- [x] `interp/triallelic.py` — `call_profile` (1/2/3 alleles, TPOX type I/II, mixture)
- [x] `interp/concordance.py` — evidence vs LongTR `discordant` flag + INEXACT propagation
- [x] `interp/profile.py` — `interpret_marker` / `interpret_run` orchestrators
- [x] CLI `frontstr interpret` end-to-end
- [x] Unit tests covering all decision branches incl. TPOX type II synthetic locus

### 1.6 HTML report — DONE

- [x] `report/payload.py` — single serializer (`serialize_run`) shared by all exports
- [x] `report/svg_charts.py` — server-side SVG: electropherogram, coverage bar, haplotype split (no JS for data)
- [x] `report/html.py` — Jinja2-based generator + SHA-256 self-stamp
- [x] `report/static/styles.css` — MinKNOW-inspired theme with dark-mode + print rules
- [x] `report/static/app.js` — vanilla JS for table sort/filter, expand-all, copy-hash
- [x] `report/templates/run_report.html.j2` — 5-section layout (Cover, Profile, QC, Loci, Audit)
- [x] CLI `frontstr report` end-to-end
- [x] Tests: payload serialization, SVG validity, full HTML parsed by lxml, deterministic (sha-stable)
- [ ] Vega-Lite + Tabulator interactive enhancement (deferred to Phase 4 polish)
- [ ] Lazy chart loading via `IntersectionObserver` (deferred — current single-file payload < 100 KB on small panels)

### 1.7 Minimal CSV/JSON exports — DONE

- [x] `exports/csv.py` — `write_profile_csv` (wide, 1×marker), `write_evidence_csv` (long, 1×cluster), `write_seqs_csv` (ISFG trail)
- [x] `exports/json.py` — `write_run_json` with `pretty` / `compact` modes
- [x] CLI `frontstr export --formats profile,evidence,seqs,json,html`
- [x] Tests: stable headers, deletion handling, tri-allelic locus serialization, JSON round-trip

**MVP Acceptance**:

1. Operator runs `frontstr run --input HG002.bam --panel codis20.yaml ...`
2. `out/S001/` contains `report.html`, `profile.csv`, `evidence.csv`, `run.json`
3. `report.html` opens in any browser, shows expected HG002 CODIS profile
4. ONT R10 simplex BAM (~1 GB) finishes in < 5 min on a laptop
5. **One real forensic lab** runs the same on their data and confirms expected genotypes.

---

## Phase 2 — Catalog, tri‑allelic polish, validation (weeks 4–6)

### 2.1 Allele catalog (3 days)

- [ ] `panel/catalog.py` — `AlleleCatalog` model
- [ ] `panel/seed_strseq.py` — STRSeq GenBank importer (PRJNA380345 + PRJNA380347)
- [ ] `interp/catalog.py` — `annotate_with_catalog` (sequence → ISFG via DB)
- [ ] `examples/catalogs/strseq_2024_06.json` — pre‑computed catalog (committable)
- [ ] Tests: known D3S1358 14a vs 14b sequences → distinct catalog hits

### 2.2 Tri‑allelic robustness (2 days)

- [ ] Per‑marker rules in panel YAML (`allow_triallelic`, `tri_balanced_thr`)
- [ ] Mixture detection: 3 peaks in non‑propensity locus → `mixture_suspected`
- [ ] Tests with TPOX type I + type II + non‑TPOX synthetic data

### 2.3 LongTR ↔ evidence concordance (1 day)

- [ ] `interp/concordance.py` — `cross_check`
- [ ] HTML report: visible discordance chip + side‑by‑side panel
- [ ] Tests: forced discordance → flag set, never silently overridden

### 2.4 Stutter overrides per marker (2 days)

- [ ] Panel YAML `stutter_overrides: {plus1A: 0.12, minus1A: 0.20}`
- [ ] Override path tested with D21S11 (compound motif)
- [ ] D3S1358 also (TCTA/TCTG mix)

### 2.5 Validation (one full week)

- [ ] HG002 + HG001 ONT R10 simplex/duplex → expected profile, F1 ≥ 0.95
- [ ] NIST 2800M control → 100% concordance with CE truth
- [ ] Mixture series 50:50, 70:30, 90:10 → mixture detected with no false `allele`
- [ ] Drop‑out study: 20×, 50×, 100×, 200× coverage → recall curves
- [ ] Write `docs/validation_report.md` with all metrics

**Phase 2 Acceptance**: Forensic lab partner runs FRONTStr on their internal
panel of known samples and signs off on a validation report.

---

## Phase 3 — Optional layers + better exports (weeks 7–8)

### 3.1 Realign to allele DB (STRspy‑style) (3 days)

- [ ] `caller/realign.py` — minimap2 against `panel_alleles.fa`
- [ ] Normalized count layer; cross‑check with evidence layer
- [ ] HTML report: extra column in evidence table

### 3.2 Phasing (2 days)

- [ ] `ingest/phase.py` — whatshap haplotag (scoped to panel ± 10 kb)
- [ ] CLI flag `--phase`
- [ ] HTML report: HP partition columns/charts active

### 3.3 Rich exports (3 days)

- [ ] `exports/xlsx.py` — multi‑sheet workbook
- [ ] `exports/vcf_extend.py` — extend LongTR VCF with `EVCOV`/`EVHP1`/…
- [ ] `exports/codis.py` — CODIS Common Message Format
- [ ] `exports/midst.py` — NIST MIDST 1.0 (mixture)
- [ ] `exports/bundle.py` — single‑click case ZIP

### 3.4 PDF (1 day)

- [ ] `report/pdf.py` — WeasyPrint render with `@media print`
- [ ] CLI flag `--export report.pdf`

---

## Phase 4 — Batch / multi‑sample (weeks 9–10)

### 4.1 Batch CLI (3 days)

- [ ] `frontstr batch --manifest batch.tsv ...`
- [ ] Parallel worker (multiprocessing) with concurrency limit
- [ ] Run‑role plumbing (`sample`, `positive_ctrl`, `negative_ctrl`, `reagent_blank`)

### 4.2 Batch report (3 days)

- [ ] `report/batch.py` — three modes: `single`, `linked`, `aggregated`
- [ ] `batch_index.html` template — heatmap, control status, alerts
- [ ] Cross‑contamination detection (shared isoalleles)
- [ ] Control validation rules (positive/negative/reagent)
- [ ] Batch ZIP layout

### 4.3 Per‑sample HTML enhancements (2 days)

- [ ] Sticky sidebar navigation (§21 of plan)
- [ ] Keyboard shortcuts
- [ ] View toggle (profile/evidence/full)
- [ ] Filter (flagged/tri/discord)

---

## Phase 5 — Web app (weeks 11–14, optional)

Only if labs need a hosted multi‑user experience.

- [ ] Postgres schema migration (`alembic`)
- [ ] FastAPI service with REST API
- [ ] React + Vite SPA (run wizard, dashboards)
- [ ] OIDC auth (Keycloak fallback to local)
- [ ] Dramatiq workers
- [ ] Docker Compose
- [ ] tus.io resumable uploads

---

## Phase 6 — Hardening + release (weeks 15–16)

- [ ] Soak test: 100 GB ONT FASTQ end‑to‑end
- [ ] Memory profiling on large batches
- [ ] Security audit (SAST + pin all deps)
- [ ] User docs, ops runbook
- [ ] `v1.0.0` tag + GitHub release + PyPI publish

---

## Dependencies snapshot (locked at Phase 1)

| Lib | Version | Why |
|---|---|---|
| `pysam` | 0.22+ | BAM I/O for pileup |
| `cyvcf2` | 0.31+ | LongTR VCF parsing |
| `edlib` | 1.3.9+ | Cluster identity |
| `pyabpoa` | 1.5+ | Cluster consensus |
| `pyyaml` | 6.0+ | Panel YAML |
| `jinja2` | 3.1+ | HTML report |
| `typer` | 0.12+ | CLI |
| `rich` | 13+ | CLI progress / logging |
| `pydantic` | 2.7+ | Panel validation |
| `pytest` | 8+ | Tests |
| `ruff` | 0.5+ | Lint |
| `mypy` | 1.10+ | Type check |

External binaries (in `PATH`):

| Binary | Version | Why |
|---|---|---|
| `minimap2` | ≥ 2.26 | Alignment |
| `samtools` | ≥ 1.18 | BAM utilities |
| `bcftools` | ≥ 1.18 | VCF utilities (Phase 4) |
| `LongTR` | ≥ 1.2 | Genotyping (Layer 1) |
| `whatshap` | ≥ 2.3 | Optional phasing |

---

## Risks & mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| `pyabpoa` install issues on macOS | Medium | High | Fall back to `spoa` subprocess; document both paths |
| LongTR VCF schema changes | Medium | High | Pin LongTR version; parse defensively; integration test on every bump |
| ONT R10 chemistry drift | Low | Medium | Quarterly validation runs; per‑chemistry profiles |
| Forensic lab rejection on validation | Medium | Critical | Engage partner early in Phase 2; iterate on feedback |
| HTML report too large on big runs | Medium | Medium | Already designed (§19 of plan): truncation + ZIP shards |
| Catalog ISFG mismatches | High | High | Manual curation step; flag low‑confidence matches |

---

## Milestones

| Milestone | When | Definition of done |
|---|---|---|
| **M0** Bootstrap | end of week 1 | CI green, `frontstr --help` works |
| **M1** MVP demo | end of week 3 | HG002 ONT BAM → HTML report with expected profile |
| **M2** Validated | end of week 6 | Forensic partner signs off on validation report |
| **M3** Feature‑complete CLI | end of week 8 | All exports + phasing + catalog working |
| **M4** Batch mode | end of week 10 | 12‑sample MinION run produces batch report |
| **M5** v1.0 GA | end of week 16 | Pinned, documented, tagged on GitHub + PyPI |
