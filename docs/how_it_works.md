# How FRONTStr works, step by step

The complete reference for what the caller does today: every stage, what it
decides, why it decides it that way, and what it does not do. Written to be
read start to finish once, then used as a lookup.

Companion documents: [`architecture.md`](architecture.md) for module layout,
[`stutter_calibration.md`](stutter_calibration.md) for the stutter measurement.

Status as of July 2026 — 402 tests passing, ~8,800 lines across 52 modules.

---

## 0. What FRONTStr is, and what it is not

**FRONTStr calls STR genotypes from long reads by itself.** Every allele, every
read count and every sequence in its output is produced by its own code, from
its own pileup of the BAM.

It is worth stating plainly because the module layout invites the opposite
conclusion. There is a `frontstr.caller` package that wraps LongTR, and it used
to be described as "Layer 1", which reads like the first step of the pipeline.
It is not. **LongTR is optional, off by default, and never contributes a
call.**

Concretely: LongTR runs only if you pass `--longtr-vcf`. When it does, the only
thing it feeds is `cross_check()`, which compares its genotype against
FRONTStr's and raises a `LONGTR_DISCORDANT` flag when they disagree. That flag
is an invitation for an analyst to look; it never changes a call. With no
`--longtr-vcf`, that step is a no-op and the result is byte-identical.

The profiles in this document were all produced on a machine with **no LongTR
installed at all**.

What FRONTStr genuinely depends on:

| Task | Provided by | When |
|---|---|---|
| Alignment (FASTQ → BAM) | `minimap2` | **Before** FRONTStr. Not part of the pipeline. |
| Reading BAM/CRAM | `pysam` | Step 1 |
| Multiple-sequence consensus | `pyspoa` (SPOA) | Step 3 |
| Edit distance | `edlib` | Step 2 |
| Everything forensic | **FRONTStr** | Steps 1–9 |
| Cross-check only | LongTR | Optional |

---

## 1. Pileup — reads to observations

**Module:** `frontstr/evidence/pileup.py` · **Input:** indexed BAM/CRAM + panel
· **Output:** one `Observation` per usable read, per marker

For each marker the panel gives a window: chromosome, start, end. The window is
the repeat array plus roughly 100 bp of flank on each side. TH01's window is
228 bp, of which the repeat array is 39 bp.

A read contributes an observation only if it **spans the entire window**, with
at least 20 bp of cleanly aligned flank on each side, and MAPQ ≥ 20. From each
qualifying read it extracts:

- the subsequence covering the window, in **reference orientation**
- the `HP` haplotype tag, when the BAM is phased
- the strand, and the mean Phred quality over the window

Two consequences worth understanding:

**This is why Illumina data does not work.** A 150 bp read cannot span a 228 bp
window. Against this panel every marker returns `no_data`. That is a design
choice, not a bug — the whole approach depends on seeing the entire locus in
one read.

**Coverage here is a count of reads, not an estimate.** Every downstream read
count traces back to a specific read that spanned the locus. Nothing is
inferred from a caller's length arithmetic.

Deletions at the window boundary are handled explicitly: a deletion covering
the start advances to the first aligned base after it, rather than dropping the
read.

---

## 2. Clustering — observations to allele candidates

**Module:** `frontstr/evidence/cluster.py`, `frontstr/motifs.py`

Two stages, and the first one is where most of the recent work went.

### 2.1 Bin by repeat-core length

The naive approach bins reads by the length of the extracted window. It fails
badly on ONT, because the window is mostly flank and ONT's dominant error is
the indel. Every read carrying a single-base indel anywhere in those ~189 bp of
flank lands in its own bin and looks like its own allele.

Measured on TH01 in HG00113, 25 reads:

```
by window length : 12 distinct values  (222, 223, 224, 227, 228, 230,
                                        235, 236, 238, 239, 240, 246)
by core length   :  4 distinct values  (28 ×9, 37 ×1, 39 ×14, 52 ×1)
```

Twelve bins for two alleles.

The core is located by `repeat_core_span`:

1. **Find every maximal run of each motif.** This catches the real array and
   also chance hits in the flanks — a stray `AATG` turns up by luck roughly
   every 256 bp.
2. **Group runs separated by ≤ 12 non-motif characters.** The threshold is not
   arbitrary. Legitimate internal spacers in forensic STRs are short: TH01's
   `ATG` is 3 bp, D21S11's structural elements (`ta`, `tca`, `tccata`) reach
   6 bp. A chance flank hit sits tens of bp away. Twelve separates the two.
3. **Keep the unit-richest group.** That is the real array; flank hits are
   discarded.

For one real TH01 read the core comes out as
`AATGAATGAATGAATGAATGAATGATGAATGAATGAATG` — note the `ATG` spacer is *inside*.
Internal structure is preserved; only flanks are trimmed.

Minus-strand markers (vWA, D5S818, CSF1PO, …) are reverse-complemented first,
because the canonical motif does not appear in reference orientation.

A read with no detectable motif run falls back to window length rather than
being dropped.

### 2.2 Why core *length* and not repeat-unit *count*

This is the decision that matters most, and the intuitive alternative is wrong:

| Allele | Structure | Repeat units | **Core length** |
|---|---|---|---|
| TH01 **9** | `[AATG]9` | 9 | **36 bp** |
| TH01 **9.3** | `[AATG]6 ATG [AATG]3` | **9** | **39 bp** |

**Both have nine units.** Binning by unit count merges them into one allele —
and 9.3 is exactly the microvariant capillary electrophoresis cannot resolve,
the thing that justifies sequencing in the first place. Core length separates
them cleanly. The same objection sinks any length-rounding scheme.

### 2.3 Merge by sequence identity within the bin

Inside each core-length bin, reads are grouped by pairwise edit-distance
identity (edlib, threshold 0.97), seed-and-grow. This second stage is what
separates iso-alleles: two alleles with the same core length but different
internal structure share a bin and are split here by sequence.

### 2.4 Result and remaining weakness

TH01 in HG00113:

```
by window : 13 clusters → [7, 5, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1]
by core   : 10 clusters → [10, 7, 1, 1, 1, 1, 1, 1, 1, 1]
```

The two real alleles went from 7 and 5 reads to **10 and 7** — matching the
final genotype `9.3(10) / 7(7)`. Across all five test samples, cluster
fragmentation dropped from **1469 to 1038 (−29%)**.

**Still not fixed:** ten clusters remain despite only four core lengths. The
14-read bin split into 10 plus four singletons. That is not the binning — it is
stage 2. The 0.97 identity threshold was set as though comparing a read to a
consensus, but `_cluster_by_identity` compares **raw read to raw read**, where
two ONT reads of the same allele diverge by ~2–4%. It does not change genotypes
here (the leftovers are singletons far below the calling threshold), but it is
the next thing to fix in this stage.

---

## 3. Consensus — clusters to sequences

**Module:** `frontstr/evidence/consensus.py`

Each cluster's member reads are collapsed into one consensus by partial-order
alignment (SPOA, global alignment). This corrects individual read errors by
voting per column.

**This step produces the single most important artefact in the pipeline.** The
ISFG string, the iso-allele match, the microvariant call and the VCF `ALT` are
all computed from this sequence. If it carries errors, everything downstream
inherits them while looking perfectly normal.

Global rather than local alignment is deliberate: clusters are length-binned,
so their members are equal-length full-window sequences, and local alignment
would be free to trim the flanks and silently change the called length — which
is the CE number.

### Why a POA backend is mandatory

Without one, the consensus degrades to *the most frequent exact sequence* —
which is a single raw read, errors included. Measured on a synthetic 252 bp
locus at ONT R10 error rates:

| Backend | 3–4 reads | 10–16 reads |
|---|---|---|
| POA | 0 edits from truth | 0 edits |
| most-common-sequence | 4–6 edits | 2–3 edits |

On the five real samples the fallback produced **four false microvariants in
202 called alleles**:

| Sample | Marker | Fallback | POA |
|---|---|---|---|
| HG00097 | TH01 | 6.3 | **7** |
| HG00154 | TH01 | 6.1 | **6** |
| HG00097 | D13S317 | 14.1 | **14** |
| HG00263 | D18S51 | 11.3 | **12** |

Every one became a clean integer. TH01 6.3 and D18S51 11.3 are not alleles that
exist in the forensic literature. In one case (HG00263 D18S51) the unpolished
consensus contained a spurious `AGAA` in the left flank — a read error that
happened to spell the motif — which corrupted the bracket structure.

`pyabpoa` does not build on macOS arm64 (it hardcodes AVX2 x86 intrinsics), so
`pyspoa` is the default. Each cluster records which method produced its
consensus; the mode fallback raises a `CONSENSUS_FALLBACK` warning on every
affected marker rather than passing silently.

---

## 4. ISFG nomenclature and allele number

**Module:** `frontstr/interp/isfg.py`, `allele_numeric.py`

From each consensus:

**The ISFG bracket string** — a greedy left-to-right scan taking the longest
motif run at each position, e.g. `[AATG]6 ATG [AATG]3`. Non-motif bases appear
in lowercase. Minus-strand markers are reverse-complemented first so the output
is in canonical orientation.

**The allele number**, by one of two routes depending on the marker:

- `period > 0` (simple markers): from length. `divmod(length − corr_value,
  period)` — integer part is full repeats, remainder is the microvariant
  decimal. TH01 39 bp → `divmod(39, 4)` = (9, 3) → **9.3**. The `corr_value` is
  the non-repeat content of the window, calibrated per marker against GRCh38.
- `period = −1` (compound markers: vWA, FGA, D21S11, …): by counting repeat
  units in the bracket string, minus a calibrated correction.

The chosen number and how it was derived both live on the model
(`Allele.number`, `Allele.number_method`), together with a single
`Allele.number_label` that every view renders. That last part matters: the CLI
used to format allele numbers independently of the report, so the same allele
could read `Δ-2` in one place and `14` in another.

---

## 5. Stutter model and classification

**Modules:** `frontstr/interp/stutter.py`, `classify.py`,
`frontstr/panel/stutter_calib.py`

Clusters above 20% of locus coverage become *parents*. For each, the model
generates the sequences its stutter products would have — the longest and
second-longest uninterrupted motif runs shortened by 1 or 2 units, or
lengthened by 1 — and an expected read count for each.

Every cluster is then classified against those expectations and against
coverage thresholds:

| Status | Meaning |
|---|---|
| `allele` | A real allele |
| `stutter` | Matches an expected stutter product at expected depth |
| `artefact` | Above the analytical threshold, below the calling threshold |
| `noise` | Below the analytical threshold (2%) |
| `inexact_allele` | Real, but the caller reconstructed rather than observed it |

### The rates are measured, not inherited

FRONTStr shipped with toaSTR-era constants: a flat 10% per LUS step, 5% per
SLUS step, forward stutter at half the reverse rate, and geometric decay for
−2. Measured against the five ONT samples (76 observations over 52 loci), all
four are wrong:

| Step | Old model | **Measured** |
|---|---|---|
| −1 | 0.100 | **0.044** |
| −2 | 0.010 (= 0.10²) | **0.011** ≈ 0.24 × the −1 rate |
| +1 | 0.050 (= 0.10 × 0.5) | **0.032** ≈ 0.73 × the −1 rate |

Forward stutter at 0.73 of reverse — rather than 0.5 — makes physical sense: on
PCR-free ONT this signal is largely sequencing indel error inside the repeat
array, which is far more symmetric than polymerase slippage.

**But the dominant effect is one the old model ignored entirely: LUS.**

| Parent LUS | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|
| Measured −1 rate | 0.010 | 0.012 | 0.035 | 0.060 | 0.122 |

More than 10× across the range where most CODIS alleles sit. No constant can
represent that; a flat 0.10 over-suppresses candidates at short-LUS loci and
under-suppresses at long-LUS ones simultaneously.

The shipped model is log-linear in LUS, R² 0.965:

```
rate(−1)   = exp(−12.1125 + 0.7159 × clamp(LUS, 10, 14))
rate(step) = rate(−1) × {−1: 1.0, −2: 0.242, +1: 0.726}
```

A plain linear fit was tried and rejected: the rates accelerate with LUS, so a
straight line fits poorly (R² 0.29) and crosses zero around LUS 9 — inside the
range where stutter is still observed, which would mean no stutter model at all
for short-LUS loci.

> **Carry this caveat.** The model is fitted on **PCR-free WGS**. It contains no
> PCR slippage component, because there was no PCR. Amplicon casework will have
> more stutter than it predicts, and it must be re-fitted before that use.
> `StutterModel.protocol` records the regime so it cannot be lost.
>
> Re-fit with: `frontstr calibrate-stutter --bam … --panel … --protocol …`

---

## 6. Haplotype-aware phantom suppression

**Module:** `frontstr/interp/haplotype.py`

The biological invariant: **a diploid locus carries exactly one allele per
haplotype.** So when two candidates are both confidently assigned to the same
haplotype and differ by less than one repeat unit, at most one can be real.

The weaker is demoted to `HP_PHANTOM`, its reads recorded on the owner's
`n_reads_absorbed`, and an `HP_PHANTOM_COLLAPSED` flag raised so the decision
is visible.

It is deliberately conservative — it fires only when every one of these holds:

- both clusters carry ≥ 3 phased reads at ≥ 80% haplotype purity (so it is a
  complete no-op on unphased BAMs)
- they share that haplotype
- they differ by ≤ 2 bp, well under one repeat unit
- the marker does not set `allow_triallelic` (a duplication genuinely puts two
  alleles on one haplotype)

### Why a read-count floor cannot do this job

The five test samples used to raise five false `mixture_suspected` flags. The
tempting fix — require a 3rd allele to clear more reads — fails on GM19038
D12S391:

```
289 bp · 10 reads · HP1        ← real allele
288 bp ·  6 reads · HP1        ← phantom of the above
301 bp ·  5 reads · HP2        ← the OTHER real allele
302 bp ·  4 reads · HP2        ← phantom of that one
```

The two strongest clusters are both phantoms of one haplotype, and the genuine
second allele is the 5-read cluster. Raising the floor deletes a real allele
and emits a confidently wrong homozygote. The haplotype tag is orthogonal
evidence that resolves it correctly, to `289 / 301`.

**With core binning in place this is now a backstop, not the workhorse:** it
fires at 1 locus out of ~110, down from 13. Core binning alone takes false
mixtures from 5 to 0.

---

## 7. Iso-allele catalog (optional)

**Module:** `frontstr/interp/catalog.py`

With `--catalog`, each called allele's repeat core is compared by edit distance
against curated sequences. An exact match adopts the published ISFG string and
the iso-allele suffix; within 2 edits it is recorded as approximate; beyond
that the live-computed values stand.

This is what turns "D3S1358 allele 14" into "**14b**" — same repeat count,
distinguishable sequence. Verified end-to-end on HG00113.

The shipped catalog (`examples/catalogs/demo_seed.json`) is a four-entry
hand-curated demonstration seed, **not** a real STRSeq import. The GenBank
fetch is still a stub.

---

## 8. Genotype call

**Module:** `frontstr/interp/triallelic.py`

Surviving candidates, sorted by read count, become a genotype:

| Candidates | Outcome |
|---|---|
| 0 | `no_data` |
| 1, or 2nd below 40% of the 1st | `homozygous` |
| 2 | `heterozygous` |
| 3+, marker allows triallelic | `triallelic_type_I` / `type_II` / `review` by balance |
| 3+, marker does not | `heterozygous`, or `mixture_suspected` if the 3rd is substantial |

A 3rd candidate must clear both the fractional calling threshold and an
absolute floor of 5 reads. With core binning and haplotype suppression both in
place, this floor is now redundant belt-and-braces.

**Amelogenin takes a separate path** (`interp/amel.py`): it is not a tandem
repeat. It counts reads at the AMELX and AMELY regions and reports X, Y or
both.

---

## 9. QC flags

**Module:** `frontstr/interp/qc.py`

Run-level conditions that depend on a laboratory threshold, applied after every
marker is called:

| Flag | Fires when | Severity |
|---|---|---|
| `DROPOUT` | No allele called at the locus | warn |
| `LOW_COVERAGE` | Called below the read floor | warn |
| `STRAND_BIAS` | An allele's strand ratio is skewed beyond chance | warn |
| `INEXACT_ALLELE` | A called allele is a reconstruction | info |
| `CE_NOMENCLATURE_OFFSET` | The number is knowingly not the kit designation | warn |

Plus the flags raised earlier in the pipeline: `MIXTURE_SUSPECTED`,
`TRIALLELIC`, `ISOALLELE`, `CONSENSUS_FALLBACK`, `HP_PHANTOM_COLLAPSED`,
`LONGTR_DISCORDANT`.

**The coverage floor is derived, not chosen.** Take the most unbalanced
heterozygote the caller still accepts (PHR 0.4, so the minor allele is 28.6% of
reads) and ask how often it falls below the calling threshold:

| Coverage | 10 | 12 | 15 | **20** | 25 | 30 |
|---|---|---|---|---|---|---|
| Dropout risk | 3.5% | 10.2% | 4.5% | **1.1%** | 1.3% | 0.3% |

20 is where the risk settles at ~1%. A first draft used 30 and flagged 12 of 25
markers on HG00113 — a flag that fires on half the panel trains reviewers to
ignore it.

**Strand bias** uses an exact two-sided binomial test written out in six lines
rather than adding SciPy. Alleles below 10 reads are not tested at all, because
the test cannot reach significance even for a perfect 0/n split, so testing
would only produce false reassurance. On 151 real called alleles it fires once
— about the rate expected by chance at p < 0.01, which says ONT R10 strand
balance is genuinely good.

**`CE_NOMENCLATURE_OFFSET` is curated in the panel**, via
`System.kit_nomenclature_note`, set for vWA and D21S11. Which markers diverge
from a kit convention is a property of the kit a laboratory compares against,
not something the code can detect.

---

## 10. Outputs and audit trail

`frontstr export --formats …` writes any of:

| Format | Contents |
|---|---|
| `profile` | CSV, one row per marker — the genotype table |
| `evidence` | CSV, one row per cluster including uncalled ones |
| `seqs` | CSV, one row per called allele with ISFG and consensus |
| `json` | The canonical record. Everything, including the audit block. |
| `html` | Self-contained offline report |
| `xlsx` | Five-sheet review workbook |
| `vcf` | Native sequence-resolved VCF (**needs `--reference`**) |

Across many samples, `frontstr tidy` flattens run JSONs into one long table
(one row per sample x marker x allele) as CSV and Parquet. `frontstr batch`
emits it automatically.

Plus `frontstr.log.jsonl`, a per-run process log.

### The VCF is native and sequence-resolved

Not an annotation of LongTR's output — a format that only exists when another
caller ran cannot serve as the benchmark interchange format.

`ALT` carries the allele's **sequence**, so iso-alleles remain distinct
records. `FORMAT/MC` carries the repeat count, `FORMAT/AD` the integer
per-allele coverage. QC warnings become `FILTER` entries; informational flags
go to `INFO/NOTE`, because resolving an iso-allele is the tool working, not a
reason to filter a call.

`--reference` is mandatory. A VCF whose REF is not the reference sequence is
not a VCF, and REF cannot be derived from a pileup.

Every per-allele `FORMAT` field is `Number=R` in REF-then-ALT order. This is
not cosmetic: emitted in call order (by depth) they misalign against `AD`
whenever a call equals the reference. Live example from TH01, where the `7`
allele *is* the reference — `MC=('9.3','7')` against `AD=(7,10)`. A benchmark
joining those columns reads every such locus wrong, and nothing about the file
looks malformed.

Verified bgzip/tabix-indexable and bcftools-queryable; REF matches GRCh38 at
all 24 emitted markers.

### Audit record

Embedded in the canonical JSON and rendered on the report's audit page: tool
version, POA backend, stutter model and protocol, every threshold that moved a
call, input hashes, the flag census, and the markers needing review.

`flags_checked` lists every code the pipeline *can* raise, so a code absent
from the counts provably means "checked and not found" rather than "never
looked at". The record is sealed with a SHA-256 over its own canonical form —
tamper evidence, not a signature.

---

## 11. Reference profile — HG00113

1000 Genomes GBR, male. ONT R10.4.1/LSK114, Dorado. This is the regression
benchmark, asserted by `tests/test_regression_hg00113.py`.

| Marker | Genotype (reads) | Coverage | Flags |
|---|---|---|---|
| D3S1358 | 14(31) | 35 | |
| vWA | 13(17) / 17(16) | 38 | `ce_nomenclature_offset` |
| FGA | 24(15) / 21(11) | 41 | |
| D8S1179 | 10(13) / 13(12) | 28 | |
| D21S11 | 35(16) / 33(10) | 29 | `ce_nomenclature_offset` |
| D18S51 | 14(11) / 15(10) | 37 | |
| D5S818 | 11(17) / 13(14) | 35 | |
| D13S317 | 11(12) / 9(11) | 28 | |
| D7S820 | 8(14) / 10(13) | 30 | |
| D16S539 | 11(20) / 8(13) | 36 | |
| TH01 | **9.3(10)** / 7(7) | 25 | |
| TPOX | 9(19) / 11(16) | 44 | |
| CSF1PO | 10(14) / 12(14) | 32 | |
| D2S1338 | 20(17) | 23 | |
| D19S433 | 14(22) / 12(17) | 45 | |
| D10S1248 | 13(19) / 15(11) | 33 | |
| D1S1656 | 13(13) / 11(12) | 27 | |
| D2S441 | 13(17) / 10(15) | 42 | |
| D12S391 | 21(20) / 20(13) | 37 | |
| D22S1045 | 16(22) | 29 | |
| AMEL | X(12) / Y(14) | 26 | |
| DYS391 | 10(14) | 16 | `low_coverage` |
| DYS393 | 14(17) | 21 | |
| DXS7132 | 15(9) | 11 | `low_coverage` |
| DXS8378 | 10(10) | 10 | `low_coverage` |

Two deviations are expected and encoded as such in the regression test:

- **vWA** reports 13/17; the kit designation is 14/16.
- **D21S11** reports 35/33; the kit designation is 29/31.

Both report the sequence-derived repeat count. For these two compound loci that
count does not equal the legacy CE-kit designation, and no single correction
reconciles them — vWA's offsets run in *opposite directions* for its two
alleles. The ISFG strings are structurally correct in both cases; only the
number differs, and both raise `CE_NOMENCLATURE_OFFSET`.

`D2S1338` reports one allele; the reference is 20/17. The second allele does
not clear the calling threshold in this slice. That is coverage, not
nomenclature.

---

## 12. Current state

### Works end to end

Pileup, clustering, POA consensus, ISFG, allele numbering, stutter,
classification, haplotype suppression, catalog annotation, genotype calling, QC
flags, all seven export formats, the audit trail, and batch mode over a
manifest. Regression-tested against a real ONT sample; CI green on ruff, ruff
format, mypy and 402 tests.

### Does not work / does not exist

- **`frontstr run` is a stub.** The FASTQ → alignment → report path is not
  wired; `ingest.align` raises. Align externally with minimap2 and start from a
  BAM.
- **No laboratory validation.** No forensic partner has signed off. No mixture
  series, no dropout study, no NIST control. The reference profile it is tested
  against comes from HipSTR on matched Illumina data — caller-vs-caller
  concordance, not validation against a reference method.
- **No CODIS `.cmf`, MIDST or PDF export.**
- **No process log in `batch` mode** — worker processes need per-process
  logging setup. The audit record *is* produced per sample.
- **Illumina data does not work** against this panel, by design.

### Known gaps in what does work

| Gap | Consequence |
|---|---|
| Identity threshold 0.97 compares raw read to raw read | Real alleles still split into singletons; no genotype effect observed, but the parameter is wrong |
| Analytical / calling thresholds (0.02 / 0.10) not data-derived | Chosen defaults. The coverage floor and stutter model *are* derived; these are not. |
| Stutter model fitted on PCR-free WGS | Under-predicts stutter on amplicon casework |
| Stutter calibrated only over LUS 10–14 | Outside that range the model clamps rather than extrapolating |
| vWA / D21S11 numbers ≠ kit designation | Cannot be compared directly against a CE profile |
| Catalog is a 4-entry demo seed | Iso-allele calling works, but has almost nothing to match against |
| X markers in males reported `homozygous` | They are hemizygous; the label is misleading |

### Priorities

1. **Widen the stutter calibration** with more R10 + Dorado samples, to extend
   the usable LUS range past 10–14.
2. **Recalibrate the identity threshold**, or switch to consensus-refined
   reassignment (seed → POA consensus → reassign reads).
3. **Derive the analytical/calling thresholds from data**, closing the last
   un-measured parameters.
4. **Either implement `frontstr run` or remove it.**
