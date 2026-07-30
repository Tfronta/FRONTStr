# How FRONTStr works, step by step

The pipeline as a sequence: each stage, what it decides, and why it decides it
that way where that is not obvious from one module. Read once, then used as a
lookup.

For running the tool, see [`../README.md`](../README.md). For the stutter
measurement, [`stutter_calibration.md`](stutter_calibration.md). For module
layout, [`architecture.md`](architecture.md).

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
repo, argv construction, BED emission, VCF parsing, all tested, because those
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

## 1. Pileup: reads to observations

**Module:** `frontstr/evidence/pileup.py` · **Input:** indexed BAM/CRAM + panel
· **Output:** one `Observation` per usable read, per marker

For each marker the panel gives a window: chromosome, start, end. The window is
the repeat array plus roughly 100 bp of flank on each side. TH01's window is
228 bp, of which the repeat array is 39 bp.

A read contributes an observation only if it **spans the entire window**, with
at least 20 bp of cleanly aligned flank on each side, and MAPQ ≥ 20. From each
qualifying read it extracts:

- the subsequence covering the window, in **reference orientation**
- the `HP` haplotype tag **and its `PS` phase block**, when the BAM is phased.
  `HP` is meaningless without the block it belongs to (§6)
- the strand, and the mean Phred quality over the window

Two consequences worth understanding:

**This is why Illumina data does not work.** A 150 bp read cannot span a 228 bp
window. Against this panel every marker returns `no_data`. That is a design
choice, not a bug, the whole approach depends on seeing the entire locus in
one read.

**Coverage here is a count of reads, not an estimate.** Every downstream read
count traces back to a specific read that spanned the locus. Nothing is
inferred from a caller's length arithmetic.

Deletions at the window boundary are handled explicitly: a deletion covering
the start advances to the first aligned base after it, rather than dropping the
read.

---

## 2. Clustering: observations to allele candidates

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
   also chance hits in the flanks, a stray `AATG` turns up by luck roughly
   every 256 bp.
2. **Group runs separated by ≤ 12 non-motif characters.** The threshold is not
   arbitrary. Legitimate internal spacers in forensic STRs are short: TH01's
   `ATG` is 3 bp, D21S11's structural elements (`ta`, `tca`, `tccata`) reach
   6 bp. A chance flank hit sits tens of bp away. Twelve separates the two.
3. **Keep the unit-richest group.** That is the real array; flank hits are
   discarded.

For one real TH01 read the core comes out as
`AATGAATGAATGAATGAATGAATGATGAATGAATGAATG`, note the `ATG` spacer is *inside*.
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

**Both have nine units.** Binning by unit count merges them into one allele,
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

The two real alleles went from 7 and 5 reads to **10 and 7**, matching the
final genotype `9.3(10) / 7(7)`. Across all five test samples, cluster
fragmentation dropped from **1469 to 1038 (−29%)**.

**Still not fixed:** ten clusters remain despite only four core lengths. The
14-read bin split into 10 plus four singletons. That is not the binning, it is
stage 2. The 0.97 identity threshold was set as though comparing a read to a
consensus, but `_cluster_by_identity` compares **raw read to raw read**, where
two ONT reads of the same allele diverge by ~2–4%. It does not change genotypes
here (the leftovers are singletons far below the calling threshold), but it is
the next thing to fix in this stage.

---

## 3. Consensus: clusters to sequences

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
would be free to trim the flanks and silently change the called length, which
is the CE number.

### Why a POA backend is mandatory

Without one, the consensus degrades to *the most frequent exact sequence*,
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
consensus contained a spurious `AGAA` in the left flank, a read error that
happened to spell the motif, which corrupted the bracket structure.

`pyabpoa` does not build on macOS arm64 (it hardcodes AVX2 x86 intrinsics), so
`pyspoa` is the default. Each cluster records which method produced its
consensus; the mode fallback raises a `CONSENSUS_FALLBACK` warning on every
affected marker rather than passing silently.

---

## 4. ISFG nomenclature and allele number

**Module:** `frontstr/interp/naming.py`, `isfg.py`, `allele_numeric.py`

From each consensus:

**The ISFG bracket string**, a greedy left-to-right scan taking the longest
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
  ranges hug the repeat array, TPOX's starts on its first base, so an aligner
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
  period)`, integer part is full repeats, remainder is the microvariant
  decimal. TH01 39 bp → `divmod(39, 4)` = (9, 3) → **9.3**. The `corr_value` is
  the non-repeat content of the window, calibrated per marker against GRCh38.
- `period = −1` (compound markers: vWA, FGA, D21S11, …): by counting repeat
  units in the bracket string, minus a calibrated correction. This is the
  arithmetic that was measurably wrong, see §12.

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
generates the sequences its stutter products would have, the longest and
second-longest uninterrupted motif runs shortened by 1 or 2 units, or
lengthened by 1, and an expected read count for each.

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

Forward stutter at 0.73 of reverse, rather than 0.5, makes physical sense: on
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
straight line fits poorly (R² 0.29) and crosses zero around LUS 9, inside the
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

It is deliberately conservative, it fires only when every one of these holds:

- both clusters carry ≥ 3 phased reads at ≥ 80% haplotype purity (so it is a
  complete no-op on unphased BAMs)
- they share that haplotype
- they differ by ≤ 2 bp, well under one repeat unit
- the marker does not set `allow_triallelic` (a duplication genuinely puts two
  alleles on one haplotype)

### Why a read-count floor cannot do this job

The five test samples used to raise five false `mixture_suspected` flags. The
tempting fix, require a 3rd allele to clear more reads, fails on GM19038
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
next are unrelated labels, the phasing tool started numbering again, not
switched chromosome. Every rule above therefore keys on `(phase_set, hp)`, and
a cluster whose tagged reads span more than one block is treated as **unphased**
rather than trusted, however pure its HP labels look.

This is not hypothetical. Across the five ONT slices, **3 of 125 loci** have
spanning, phased reads from more than one block. HG00097 D13S317 is the clear
case: a 14-read cluster reporting 100% HP2, drawn 4 reads from one block and 3
from another.

No call changes: those loci are balanced enough that no haplotype rule fires,
so this is a latent-evidence fix, not a genotype fix. The affected loci raise
`PHASE_BLOCK_SPLIT` (WARN) and the trace marks the cluster
`[2 phase blocks, HP not comparable]`, because a reviewer comparing HP counts
by eye would otherwise draw the conclusion the caller deliberately refused to.

A BAM with `HP` but no `PS` keeps the pre-`PS` behaviour: refusing there would
silently disable haplotype reasoning for every phasing tool that omits the tag.

### The same invariant, read backwards

If two candidates on the *same* haplotype cannot both be real, then two on
*opposite* haplotypes cannot be the same allele, however unbalanced their read
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

This is what turns "D3S1358 allele 14" into "**14b**", same repeat count,
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
consensus of `17/20`, a false exclusion, the costliest error this caller can
make. At ONT depths a 5-vs-17 split is ordinary sampling; the phasing says
plainly that these are two alleles.

So the ratio now yields to `on_opposite_haplotypes()` (§6), under the same gates
 ≥ 3 tagged reads at ≥ 80% purity on *both* sides, hence a no-op on unphased
BAMs, where the ratio is again the only evidence there is. The rescued allele
carries `hp_rescued` and the marker raises `HP_RESCUED_HET`, so a call resting
on phasing rather than balance is never silent.

It stays rare by construction: across the five test samples it fires at exactly
one locus. A version of this that fires broadly is no longer an appeal to
phasing evidence, it is a lowered threshold, `test_phasing_rescue_does_not_invent_heterozygotes`
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
markers on HG00113, a flag that fires on half the panel trains reviewers to
ignore it.

**Coverage is the evidence behind the call.** `Cov` in the CLI and the report
is `called_reads`, the reads supporting the reported genotype, not every read
spanning the window. Reads the caller rejected were rejected for a reason;
counting them overstates the evidence and makes the per-allele numbers look as
if they fail to add up. `total_reads` survives internally because it is the
denominator every fraction threshold is measured against, and any spanning read
that supported no called allele shows as a trailing `+n`.

`low_coverage` moved to the same number, and that is a correctness fix rather
than tidying: the binomial models *N reads split between two alleles*, and
reads clustering into neither are not draws from that pair. Measuring the
spanning total instead let HG00263 D18S51 pass silently, called on 11 reads
out of 33 spanning, missing the second allele Illumina sees.

Re-deriving the floor against ONT depths confirmed **20** but corrected the
claim attached to it. Read as "from this coverage upward the risk never again
exceeds", the curve is a sawtooth, and the trade is:

| Floor | Dropout risk | Loci flagged |
|---|---|---|
| 17 | 9.7% | 12% |
| **20** | **5.7%** | **25%** |
| 25 | 4.6% | 43% |

20 is the knee: 17 → 20 halves the risk for 13 points of flag rate, while
20 → 25 buys one point of risk for eighteen. The docstring previously claimed
~1.1% risk at N=20; that figure came from measuring against the spanning total
and does not survive. A quarter of loci flagged at ~30x ONT is the expected
rate, not a symptom, median called coverage in the reference slices is 26 and
the lower quartile is 20, so the floor sits at Q1 by construction.

**Strand bias** uses an exact two-sided binomial test written out in six lines
rather than adding SciPy. Alleles below 10 reads are not tested at all, because
the test cannot reach significance even for a perfect 0/n split, so testing
would only produce false reassurance. On 151 real called alleles it fires once
 about the rate expected by chance at p < 0.01, which says ONT R10 strand
balance is genuinely good.

**`CE_NOMENCLATURE_OFFSET` is curated in the panel**, via
`System.kit_nomenclature_note`, set for vWA and D21S11. Which markers diverge
from a kit convention is a property of the kit a laboratory compares against,
not something the code can detect. The flag is now suppressed when every called
allele was named by STRNaming, because in that case the number *is* the kit
designation, warning a reviewer off comparing it would cause exactly the false
exclusion the flag exists to prevent. The notes stay in the panel because they
still describe the fallback path.

---

## 10. Where the rest is

The three sections that used to close this document were re-explaining things
that are now documented where they belong, and two copies of an explanation
drift apart. What they covered, and where it lives:

| Was here | Is now in |
|---|---|
| `--log`, `--trace`, the profile table, export formats, the VCF, `--bed`, `doctor` | [`../README.md`](../README.md), which is where someone looks before running anything |
| The parameter table and what overriding one costs | `frontstr/params.py`, beside each default, and summarised in the README |
| The coverage floor and balance band derivations | `frontstr/interp/qc.py`, beside the constants they justify |
| The HG00113 reference profile | [`../demodata/README.md`](../demodata/README.md), beside the sample it describes |
| Current state and known gaps | the README's Limitations, and [`../ROADMAP.md`](../ROADMAP.md) |
| The audit trail's four layers | `frontstr/audit.py` |

This file keeps what none of those can carry: the pipeline as a sequence, with
the reasoning for each stage's design where that reasoning is not obvious from
one module.
