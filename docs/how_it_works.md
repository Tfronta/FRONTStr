# How FRONTStr works, step by step

The complete reference for what the caller does today: every stage, what it
decides, why it decides it that way, and what it does not do. Written to be
read start to finish once, then used as a lookup.

Companion documents: [`architecture.md`](architecture.md) for module layout,
[`stutter_calibration.md`](stutter_calibration.md) for the stutter measurement.

Status as of July 2026 — 591 tests passing, ~11,200 lines across 56 modules.

---

## 0. What FRONTStr is, and what it is not

**FRONTStr calls STR genotypes from long reads by itself.** Every allele, every
read count and every sequence in its output is produced by its own code, from
its own pileup of the BAM.

It is worth stating plainly because the module layout invites the opposite
conclusion. There is a `frontstr.caller` package that wraps LongTR, and it used
to be described as "Layer 1", which reads like the first step of the pipeline.
It never was, and as of July 2026 it is not connected at all: **no external
caller runs, and no command exposes one.**

It used to be an optional cross-check behind `--longtr-vcf` that raised a flag
on disagreement without ever changing a call. That was unwired, along with the
`interp.concordance` glue it fed and the `call` command that ran it. The
reasoning is that benchmarking FRONTStr against another caller is a different
activity from running FRONTStr, and doing both in one command made it hard to
say which tool produced a number. The `frontstr.caller` package stays in the
repo — argv construction, BED emission, VCF parsing, all tested — because those
are the tedious parts of a future benchmark harness.

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
- the `HP` haplotype tag **and its `PS` phase block**, when the BAM is phased —
  `HP` is meaningless without the block it belongs to (§6)
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

**Module:** `frontstr/interp/naming.py`, `isfg.py`, `allele_numeric.py`

From each consensus:

**The ISFG bracket string** — a greedy left-to-right scan taking the longest
motif run at each position, e.g. `[AATG]6 ATG [AATG]3`. Non-motif bases appear
in lowercase. Minus-strand markers are reverse-complemented first so the output
is in canonical orientation.

**The allele number**, from STRNaming whenever it is available. The ISFG DNA
Commission (Gettings et al. 2024) names STRNaming as *the program* that produces
the designation, rather than a rule set to reimplement, and FRONTStr follows
that: the number is whatever STRNaming reports for the marker's standard
reporting range (`CE29_TCTA[4]TCTG[6]…` → 29).

Two details make this work on long reads:

- The reporting range is located in the consensus by **aligning the reference
  flanks that sit just outside it**, not by reference coordinates. STRNaming's
  ranges hug the repeat array — TPOX's starts on its first base — so an aligner
  that places a long allele's extra units outside the boundary would otherwise
  make them vanish. Coordinate slicing named both HG00113 TPOX alleles CE8 when
  the truth is 9/11.
- Reference sequence comes from a **committed GRCh38 slice cache**
  (`frontstr/interp/data/strnaming_ranges.tsv`, ~13 kB, built by
  `frontstr/panel/seed_strnaming.py`). No network call and no `--reference`
  requirement, so naming is reproducible byte-for-byte.

The legacy arithmetic remains as the fallback for markers STRNaming defines no
range for (DYS393, AMEL), and for any consensus the range cannot be located in:

- `period > 0` (simple markers): from length. `divmod(length − corr_value,
  period)` — integer part is full repeats, remainder is the microvariant
  decimal. TH01 39 bp → `divmod(39, 4)` = (9, 3) → **9.3**. The `corr_value` is
  the non-repeat content of the window, calibrated per marker against GRCh38.
- `period = −1` (compound markers: vWA, FGA, D21S11, …): by counting repeat
  units in the bracket string, minus a calibrated correction. This is the
  arithmetic that was measurably wrong — see §12.

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

### Haplotype labels are local to a phase block

`HP` only means something inside its own `PS`. HP1 in one block and HP1 in the
next are unrelated labels — the phasing tool started numbering again, not
switched chromosome. Every rule above therefore keys on `(phase_set, hp)`, and
a cluster whose tagged reads span more than one block is treated as **unphased**
rather than trusted, however pure its HP labels look.

This is not hypothetical. Across the five ONT slices, **3 of 125 loci** have
spanning, phased reads from more than one block. HG00097 D13S317 is the clear
case: a 14-read cluster reporting 100% HP2, drawn 4 reads from one block and 3
from another.

No call changes — those loci are balanced enough that no haplotype rule fires —
so this is a latent-evidence fix, not a genotype fix. The affected loci raise
`PHASE_BLOCK_SPLIT` (WARN) and the trace marks the cluster
`[2 phase blocks — HP not comparable]`, because a reviewer comparing HP counts
by eye would otherwise draw the conclusion the caller deliberately refused to.

A BAM with `HP` but no `PS` keeps the pre-`PS` behaviour: refusing there would
silently disable haplotype reasoning for every phasing tool that omits the tag.

### The same invariant, read backwards

If two candidates on the *same* haplotype cannot both be real, then two on
*opposite* haplotypes cannot be the same allele — however unbalanced their read
counts. `on_opposite_haplotypes()` exposes that reading, and §8 uses it to stop
the peak-height ratio from collapsing a genuine heterozygote. Same evidence,
same gates, opposite direction: one deletes alleles, the other rescues them.

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
| 1, or 2nd below 40% of the 1st **and not on the opposite haplotype** | `homozygous` |
| 2 | `heterozygous` |
| 3+, marker allows triallelic | `triallelic_type_I` / `type_II` / `review` by balance |
| 3+, marker does not | `heterozygous`, or `mixture_suspected` if the 3rd is substantial |

A 3rd candidate must clear both the fractional calling threshold and an
absolute floor of 5 reads. With core binning and haplotype suppression both in
place, this floor is now redundant belt-and-braces.

### The 40% floor is not the best evidence available

`min_phr_for_het` is inherited from capillary electrophoresis, where peak height
is all there is. On phased long reads it is not. HG00113 D2S1338:

```
20 · 17 reads · 100% HP1
17 ·  5 reads · 100% HP2     PHR = 0.29, under the 0.4 floor
```

That was called homozygous `20`, against an Illumina, LongTR **and** STRspy
consensus of `17/20` — a false exclusion, the costliest error this caller can
make. At ONT depths a 5-vs-17 split is ordinary sampling; the phasing says
plainly that these are two alleles.

So the ratio now yields to `on_opposite_haplotypes()` (§6), under the same gates
— ≥ 3 tagged reads at ≥ 80% purity on *both* sides, hence a no-op on unphased
BAMs, where the ratio is again the only evidence there is. The rescued allele
carries `hp_rescued` and the marker raises `HP_RESCUED_HET`, so a call resting
on phasing rather than balance is never silent.

It stays rare by construction: across the five test samples it fires at exactly
one locus. A version of this that fires broadly is no longer an appeal to
phasing evidence, it is a lowered threshold — `test_phasing_rescue_does_not_invent_heterozygotes`
guards that.

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
not something the code can detect. The flag is now suppressed when every called
allele was named by STRNaming, because in that case the number *is* the kit
designation — warning a reviewer off comparing it would cause exactly the false
exclusion the flag exists to prevent. The notes stay in the panel because they
still describe the fallback path.

---

## 10. Outputs and audit trail

### Watching a run, not just reading its conclusion

`frontstr interpret --log` prints the process log to **stderr** as it goes: the
run configuration, then one line per marker with the call rule, the coverage,
the cluster count and how each allele number was derived.

```
[info ] run.start      bam=HG00113.codis.bam n_markers=25 min_mapq=20 strnaming=True
[debug] marker.called  marker=vWA call_rule=heterozygous alleles=['14', '16']
                       n_clusters=7 total_reads=38 number_method=['strnaming']
```

It goes to stderr so `frontstr interpret … --log > table.txt` still separates
the result table from the commentary. The same events land in
`frontstr.log.jsonl` as JSON when `export` writes an output directory — the two
sinks share one event stream and render it differently, because a file wants to
be diffed and a terminal wants to be read.

`naming.fallback` lines appear for clusters STRNaming declined to name:

```
[debug] naming.fallback  marker=FGA reason=anchor_low_identity n_reads=1
                         consensus_method=single number_method=bracket_count
```

These are **not errors**. Every one on the reference slices is a one-read
cluster whose unpolished ONT indels mangled a 30 bp flank past the identity
guard — the guard refusing to name a corrupt sequence, which is what it exists
for. Called alleles use a POA consensus and are all named. The read count and
`consensus_method=single` are in the line precisely so this is judgeable at a
glance rather than alarming.

The same condition on a **called** allele is logged as
`naming.fallback_on_called_allele` at *warning* level, because then a reported
number came from the legacy arithmetic that is wrong at six markers.

### `--trace`: the whole locus, step by step

`--log` says what was concluded. `--trace` says how. It narrates every stage in
pipeline order, for every locus, so a genotype can be followed back to the
reads — which is what makes a call reviewable rather than merely reported.

```
── FRONTStr 0.1.0.dev0
  Detected                                1 BAM
      tests/data/ont_slices/HG00113.codis.bam
  Panel                                   CODIS 20 + sex 2026.06-ont-wide-sex
  Markers in panel                        25
  Read filters                            MAPQ >= 20, 20 bp clean flank each side
  Thresholds                              analytical 2%, calling 10%, cluster identity 0.97
  Consensus backend                       poa_spoa
  Allele naming                           STRNaming, offline slice cache, 23 markers

── TH01  chr11:2170988-2171215  (228 bp window, motif AATG, period 4)
  Reads fetched around the window         25
  Spanning the whole window               25   (total locus coverage)

  Step 1 — grouped by length              4 bins, using repeat-core length (flank indel errors ignored)
      39 bp core                          14 reads
      28 bp core                          9 reads
      37 bp core                          1 read
      52 bp core                          1 read
  Step 2 — split by sequence              10 clusters, 6 more than bins (identity below 0.97 separates)
  Step 3 — consensus per cluster          poa_spoa

  Candidates, strongest first             reads shown are per-allele coverage
    * #0    10 reads  40.0%   239 bp  HP1 0 / HP2 10 / untagged 0
        name      CE9.3_TGAA[6]TGA[1]TGAA[3]
        number    9.3   via STRNaming
        verdict   allele
    * #1     7 reads  28.0%   228 bp  HP1 7 / HP2 0 / untagged 0
        name      CE7_TGAA[7]
        number    7   via STRNaming
        verdict   allele
      #2      1 read   4.0%   235 bp  HP1 0 / HP2 1 / untagged 0
        name      CE8.2_TGAA[6]TGA[1]TGAA[3]_-1G>-_-5CA>AC_…
        consensus single  ← not polished by POA
        verdict   artefact  (above analytical, below the 10% calling threshold)
  …
  Sequences (flank … repeat core … flank):
    * #0  …GACTCCATGGTG  AATGAATGAATGAATGAATGAATGATGAATGAATGAATG   AGGGAAATAAGG…
    * #1  …GACTCCATGGTG  AATGAATGAATGAATGAATGAATGAATG              AGGGAAATAAGG…
      #2  …ACACGCACTGTG  AATGAATGAATGAATGAATGAATGATGAATGAATGAATG   GAGGGAGAGGGA…
  …
  Genotype                                9.3 (10 reads), 7 (7 reads)   [heterozygous]
  Coverage                                25 at the locus; 17 on called allele(s);
                                          8 on discarded candidates; phased
                                          9.3: HP1 0 HP2 10 / 7: HP1 7 HP2 0
```

That single block answers the questions the profile table cannot: TH01's two
core-length bins (39 bp and 28 bp) *are* the 9.3 and the 7, separated cleanly
before any clustering; the two called alleles sit on opposite haplotypes; and
the eight leftover candidates are single unpolished reads rejected by a named
threshold rather than silently dropped.

**Every discarded candidate says which rule discarded it** — `noise` cites the
analytical threshold, `artefact` the calling threshold, `stutter` the expected
count from its parent, `hp_phantom` the one-allele-per-haplotype invariant. So
does a rescue: an allele called on phasing against the peak-height ratio says
so on its verdict line.

**The three clustering steps are named, not implied.** Step 1 groups reads by
repeat-core length; step 2 splits each of those groups by sequence identity, so
two alleles of the same length but different internal structure — iso-alleles —
come apart; step 3 collapses each surviving group into one consensus. When step
2 reports the same number of clusters as step 1 had bins, no bin contained two
distinct sequences. When it reports more, the extra ones are almost always
single noisy reads that failed to join their own allele's seed.

**Per-allele coverage is on the conclusion, not just in the table.** Integer
per-allele read counts are the headline claim of this caller over a
length-based one, so the genotype line carries them, and a `Coverage` line
splits the locus total into what supports the called alleles, what went to
discarded candidates, and how each called allele divides across haplotypes.
They are also in `profile.csv` (`alleleN_cov`, `alleleN_hp1`, `alleleN_hp2`),
in `evidence.csv` one row per cluster, in the XLSX, and in the VCF `FORMAT/AD`.

**The bases are shown, aligned.** The repeat core prints in full with its
flanks cut to 12 bp on each side — a panel window is ~80% flank, and printing
all of it would bury the part that distinguishes two alleles. Cores are padded
into one column so reading downward makes a 4 bp step, an interrupted motif or
a single substitution visible without counting characters. Above, `#0`'s `ATG`
interruption against `#1`'s clean `AATG` run *is* 9.3 versus 7, and `#2`'s
different left flank marks it as a misaligned read rather than a third allele.

**A called allele's core is never truncated** — D21S11's runs past 180 bp and
still prints whole, because the point of showing bases is to validate the
genotype. Uncalled candidates are cut at 96 bp, and the cut is disclosed rather
than silent. Sequences are in canonical (motif) orientation, so a minus-strand
marker reads the same way as its ISFG string.

**Allele balance, and the QC verdict.** The genotype line carries the flags
that actually fired, and a heterozygote gets an `Allele balance` line:

```
  Genotype        20 (17 reads), 17 (5 reads)   [heterozygous]   INFO: hp_rescued_het
  Allele balance  0.77  (uneven; 0.50 is even, balanced up to 0.65, …)
```

AB is the **strongest called allele over the called pair**, so it runs from
0.50 (even) to 1.0. Stating that convention matters: with the strongest allele
on top the scale is one-sided, and a band written as if it were symmetric
around 0.5 would have half of itself unreachable. It replaces the peak-height
ratio in output rather than joining it — `AB = 1/(1+PHR)`, the same
measurement, and two numbers for one quantity is how the same allele once read
`Δ-2` in the CLI and `14` in the report.

The landmarks line up: 0.50 even, **0.65** the balanced band, **0.714** is
`min_phr_for_het` = 0.4 below which no het is called on read counts alone, and
0.773 is HG00113 D2S1338 — called het on phasing alone. So the band sits inside
the callable range and warns rather than rejects. It fires at 11 of ~200 called
loci across the five slices, all of them correctly genotyped: the flag means
*watch this locus for dropout elsewhere*, not *this call is wrong*. It is
suppressed when `HP_RESCUED_HET` already fired, since that flag says the ratio
was the problem.

**There is no aggregated `PASS`.** A single green label standing for several
checks teaches a reviewer to stop reading the individual ones, and one that
shows on almost every locus carries no information. A clean locus says nothing;
a flagged one names what fired, worst severity first. The same applies to the
`interpret` table, which gained `AB` and `QC` columns.

**The read funnel closes.** `kept + rejected == fetched` at every locus, with
each rejection under a named reason (`not a primary alignment`, `MAPQ below
threshold`, `does not reach the left flank anchor`, …). A coverage number a
reviewer cannot reconcile against the BAM is not an auditable one, so this is
asserted end to end in `test_trace_accounts_for_every_read_at_every_locus`.

A run summary closes the trace with the totals, any locus without a genotype,
and any allele that fell back to the legacy CE path.

`--trace-out FILE` writes the narrative to a file instead of stderr.

Still missing, and worth adding: the POA step reports only which backend ran,
not how many positions it corrected.

### The profile table

One row per locus, answering the whole question without a second view: sample,
both allele numbers, **per-allele coverage**, total, allele balance, the QC
flags that fired, and **both consensus sequences**.

Each sequence sits in its own fixed-width cell that scrolls sideways. A ~250 bp
consensus wrapped over ten lines turns every row into a paragraph, and
truncating it would hide the one thing a sequence-resolved caller exists to
show. The table itself scrolls as a unit so the sequence columns cannot squeeze
the rest into slivers.

Flags used to live **only** inside the expandable per-locus cards at the bottom
of the report, so a reviewer scanning the profile table could not see that a
locus was flagged at all — the XLSX export had been doing this correctly
(tinted rows plus a QC sheet) while the HTML had not. Now the row is tinted by
worst severity, the QC cell names each flag, and the toolbar filters to flagged
or clean. Typing a flag code into the search box works too.

There is no aggregated `PASS`, for the reason given above: a label that appears
on almost every row stops being read.

**One repeat string everywhere.** `Allele.repeat_label` is STRNaming's name when
the marker has a reporting range, and the legacy `compress_isfg` scan otherwise.
Every view renders it — the same reason `number_label` exists. Before this, the
trace printed `CE9.3_TGAA[6]TGA[1]TGAA[3]` while the HTML, CSV and XLSX printed
the full-window scan: a hundred lowercase flank bases before the brackets even
start. `isfg_source` (`strnaming` | `bracket_scan`) travels with it, and the raw
window scan is still available as `isfg_window` in the JSON.

### Parameters, and what changing one costs

**Module:** `frontstr/params.py`

Every knob a run can turn lives in one table with its default and where that
default came from. A run prints the whole table before it starts — not only what
was typed, because the values that decide a call are usually the ones nobody
typed, and a run listing only its command line reads as if it were barely
configured.

Provenance is the reason the table exists rather than a plain dict:

| | Meaning | Changing it |
|---|---|---|
| `derived` | Computed from measured data — `low_coverage_reads` from a binomial dropout calculation, `min_reads_third` from the known-bug #6 phantom work | **Marks the run** |
| `chosen` | Defensible but picked; the analytical and calling thresholds are the honest examples | Ordinary tuning |
| `convention` | Inherited from forensic practice — `min_phr_for_het` is the CE heterozygote-balance rule, which is why phasing overrides it | Ordinary tuning |

Overriding a `derived` default raises `NON_DEFAULT_THRESHOLD` on **every
marker**, carrying the value, the default, and why that default existed. Marker
level rather than run level because that is where flags already live and where
every consumer already looks: the audit census, the XLSX QC sheet, the HTML row
tint. A run-level-only note is the one thing a reviewer scanning per-locus rows
never sees.

The point is not to stop anyone testing a threshold — `--min-reads-third 2` is a
legitimate experiment. It is that six months later the profile still says which
numbers were changed, and that the ones with measured backing are called out
rather than buried among forty.

`tests/test_params.py` asserts the table's defaults equal the ones the pipeline
actually uses; a table that drifts is worse than none, because the report would
confidently print a default nothing is using.

### Bringing your own regions

`--bed FILE` replaces `--panel`, which is the escape hatch HipSTR and LongTR
have and a YAML-only caller did not: point the tool at your own intervals
without curating a panel first.

```
chrom  start  end  MOTIF  name
```

**The motif is required**, and that is not a parsing convenience. Without it
reads cannot be binned by repeat-core length, and binning on raw window length
instead took TH01 from 2 bins to 12. Better to refuse the file than to call from
it badly.

**Coordinates have no default anywhere in `panel/bed.py`.** Standard BED is
0-based half-open; the panel YAML and HipSTR's region files are 1-based
inclusive. One base at a window edge produces a plausible wrong answer rather
than a crash, so `--bed-coords` is explicit (`bed0`, the default for reading,
or `panel1`) and every writer states which it emits. The report's BED block and
`--bed` use the same convention, so the round trip is exact — asserted in
`tests/test_panel_bed.py`.

What a BED cannot carry is calibration: no `period`, no `corr_value`, no
`marker_type`. For markers STRNaming has a reporting range for, that costs
nothing — which is only true since STRNaming became the canonical namer, and is
what makes this feature worth having now. For the rest, the number is an
uncalibrated repeat count rather than a kit allele, so those markers get a
`kit_nomenclature_note` and every call at them raises
`CE_NOMENCLATURE_OFFSET` — the same machinery a curated panel uses to say the
same thing. The run also prints them by name before it starts.

Measured on HG00113, feeding back the panel's own windows as BED: 23 of 25
genotypes identical to the YAML run. The two that differ are exactly the two
the warning names — DYS393 (no STRNaming range) and AMEL (a BED cannot say
`marker_type: amel`, so sex typing degrades to the STR path).

### Checking the installation

`frontstr doctor` with no arguments checks the install alone: the POA backend,
the STRNaming slice cache, the compiled dependencies. Worth running after any
install, because the failure mode is quiet — without a POA backend FRONTStr
still produces a complete profile, built from unpolished single reads, and the
damage surfaces as microvariants that are not in the sample (4 in 202 called
alleles, measured). It exits non-zero when something is broken, so it can gate
a batch.

Given `--bam` and `--panel` it adds the per-marker table and a phasing line:
what fraction of reads carry `HP`, whether `PS` is present, and whether any
sampled locus spans more than one phase block. Both haplotype rules are
silently disabled on an unphased BAM, and it is better to know that before
reading a profile than to wonder why a rescue never fired.

### Provenance sections

The Raw / Audit page carries three things a reviewer would otherwise have to
take on trust:

- **Command** — the exact `sys.argv`.
- **Parameters in force** — every value the run used, defaults included. The
  command line alone is misleading, because the parameters that decide a call
  are usually the ones nobody typed.
- **Panel windows (BED)** — the extraction intervals verbatim, 1-based
  inclusive, `chrom start end motif name`, with a copy button. Naming the loci
  is not the same as showing the intervals they came from; this pastes straight
  into samtools or IGV to look at the same reads the caller saw.

All three are omitted rather than faked when a library caller supplies none of
them.

### Export formats

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

Not an annotation of another caller's output — a format that only exists when another
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
| vWA | 14(17) / 16(16) | 38 | |
| FGA | 24(15) / 21(11) | 41 | |
| D8S1179 | 10(13) / 13(12) | 28 | |
| D21S11 | 31(16) / 29(10) | 29 | |
| D18S51 | 14(11) / 15(10) | 37 | |
| D5S818 | 11(17) / 13(14) | 35 | |
| D13S317 | 11(12) / 9(11) | 28 | |
| D7S820 | 8(14) / 10(13) | 30 | |
| D16S539 | 11(20) / 8(13) | 36 | |
| TH01 | **9.3(10)** / 7(7) | 25 | |
| TPOX | 9(19) / 11(16) | 44 | |
| CSF1PO | 10(14) / 12(14) | 32 | |
| D2S1338 | 20(17) / 17(5) | 23 | `hp_rescued_het` |
| D19S433 | 14(22) / 12(17) | 45 | |
| D10S1248 | 13(19) / 15(11) | 33 | |
| D1S1656 | 13(13) / 11(12) | 27 | |
| D2S441 | 13(17) / 10(15) | 42 | |
| D12S391 | 21(20) / 20(13) | 37 | |
| D22S1045 | 16(22) | 29 | |
| AMEL | X(12) / Y(14) | 26 | |
| DYS391 | 10(14) | 16 | `low_coverage` |
| DYS393 | 14(17) | 21 | |
| DXS7132 | 14(9) | 11 | `low_coverage` |
| DXS8378 | 10(10) | 10 | `low_coverage` |

**Every marker with a Book1 reference now matches it**, vWA and D21S11
included. Both used to deviate — 13/17 and 35/33 against kit values of 14/16
and 29/31 — because the bracket count is not the kit designation and no single
`corr_value` reconciles it; vWA's offsets run in *opposite directions* for its
two alleles. Adopting STRNaming (§4) resolved both.

That deviation turned out not to be limited to those two. Naming the GRCh38
reference window itself and comparing against the published `ref_ce` shows the
old arithmetic was also off at **D2S1338, D1S1656, D2S441 and DXS7132** — for
the compound markers by an allele-structure dependent amount, which is why it
looked correct on some samples and not others. `DXS7132` moved 15 → 14 here;
its `corr_value` had been calibrated against FRONTStr's own output rather than
an external truth.

`D2S1338` is now 20/17, matching the reference. Its minor allele is carried by
5 reads against 17 — under the 40% peak-height floor — and is called because
phasing puts the two clusters on opposite haplotypes (§8).

### Concordance against the external reference

Measured across all five slices against the `Illumina` sheet of
`1000GEN-ONT-Merged-Compar.xlsx`, restricted to markers that sheet actually
calls: **91 of 95 genotypes (95.8%)**. The four remaining are all the minor
allele of a heterozygote lost to a threshold, none to nomenclature:

| Sample | Marker | Reference | FRONTStr | Why |
|---|---|---|---|---|
| HG00154 | D18S51 | 13/14 | 14 | major cluster is 11 HP1 / 7 HP2, so phasing cannot rescue it |
| HG00263 | D18S51 | 12/16 | 12 | minor allele is 3 reads = 9.1%, just under the 10% calling threshold |
| HG00154 | D5S818 | 12/13 | 13 | 11 reads at the locus; genuine low-coverage dropout |
| HG00154 | FGA | 22/19 | 24/19 | **not a FRONTStr error** — LongTR and STRspy both call 24 |

The FGA row is worth reading carefully: Illumina is the outlier there, out-voted
3–1 by the long-read methods. Short reads under-call long FGA alleles, which is
the failure mode this project exists to avoid.

---

## 12. Current state

### Works end to end

Pileup, clustering, POA consensus, ISFG, allele numbering, stutter,
classification, haplotype suppression, catalog annotation, genotype calling, QC
flags, all seven export formats, the audit trail, and batch mode over a
manifest. Regression-tested against a real ONT sample; CI green on ruff, ruff
format, mypy and 591 tests.

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
| DYS393 and AMEL have no STRNaming range | Both keep the legacy arithmetic; DYS393's number is unverified against an external truth |
| The sex markers have no external reference at all | DYS391 is `NA` in the reference workbook and DYS393/DXS7132/DXS8378 are absent; those four expectations are FRONTStr's own output |
| D21S11 has no orthogonal reference | HipSTR does not call it; 29/31 rests on LongTR and STRspy, both long-read |
| Two heterozygotes still lost to thresholds | HG00154 D18S51 (phasing too impure to rescue) and HG00263 D18S51 (3 reads, just under the 10% calling threshold) |
| The bracket string is still FRONTStr's own, not STRNaming's | ISFG Recommendation 2 asks for STRNaming's formatting; only the *number* has been adopted so far |
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
