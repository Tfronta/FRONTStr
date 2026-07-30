# How FRONTStr works, step by step

The pipeline as a sequence: each stage, what it decides, and why it decides it
that way where that is not obvious from one module. Read once, then used as a
lookup.

For running the tool, see [`../README.md`](../README.md). For the stutter
measurement, [`stutter_calibration.md`](stutter_calibration.md). For module
layout, [`architecture.md`](architecture.md).

---

## 0. Scope

Input: an indexed BAM or CRAM aligned to GRCh38, and a panel.
Output: one genotype per marker, with per-allele read depth, ISFG nomenclature
and QC flags.

Genotypes are called from FRONTStr's own pileup of the BAM. No external caller
runs and no command exposes one.

| Task | Provided by | When |
|---|---|---|
| Alignment (FASTQ to BAM) | `minimap2` | Before FRONTStr, not part of the pipeline |
| Reading BAM/CRAM | `pysam` | Stage 1 |
| Edit distance | `edlib` | Stage 2 |
| Multiple-sequence consensus | `pyspoa` | Stage 3 |
| Allele nomenclature | `strnaming` | Stage 4 |

---

## 1. Pileup: reads to observations

**Module:** `frontstr/evidence/pileup.py`
**Input:** indexed BAM/CRAM and a panel. **Output:** one `Observation` per
usable read, per marker.

Each marker defines a window: chromosome, start, end. The window is the repeat
array plus about 100 bp of flank each side. TH01's window is 228 bp, of which
39 bp is the array.

A read contributes an observation only if it spans the entire window, with at
least `--flank-anchor` bp of cleanly aligned flank on each side and MAPQ at or
above `--min-mapq`. From each qualifying read FRONTStr extracts:

- the subsequence covering the window, in reference orientation
- the `HP` haplotype tag and its `PS` phase block, when the BAM is phased.
  `HP` is meaningless without the block it belongs to (stage 6)
- the strand, and the mean Phred quality over the window

Read counts downstream are counts of these observations. Nothing is inferred
from length arithmetic.

Deletions at a window boundary advance to the first aligned base after the
deletion rather than dropping the read.

**Illumina data does not work against this panel.** A 150 bp read cannot span a
228 bp window, so every marker returns `no_data`. The approach requires the
whole locus in one read.

---

## 2. Clustering: observations to allele candidates

**Module:** `frontstr/evidence/cluster.py`, `frontstr/motifs.py`

Reads are grouped in two stages: by repeat-core length, then by sequence within
each group.

### 2.1 Bin by repeat-core length

Reads are binned on the length of the repeat core, not of the extracted window.

A panel window is about 80% flank, and ONT's dominant error is the indel. On
window length, one stray flank indel puts a read in its own bin. TH01 in the
demo sample, 25 reads:

| Binned on | Distinct bins |
|---|---|
| Window length | 12 (222, 223, 224, 227, 228, 230, 235, 236, 238, 239, 240, 246) |
| Core length | 4 (28 ×9, 37 ×1, 39 ×14, 52 ×1) |

`repeat_core_span` locates the core:

1. Find every maximal run of each motif. This catches the array and also chance
   hits in the flanks, where a stray `AATG` occurs roughly every 256 bp.
2. Group runs separated by 12 or fewer non-motif characters. Internal spacers
   in forensic STRs are short: TH01's `ATG` is 3 bp, D21S11's `ta`, `tca` and
   `tccata` reach 6 bp. Chance flank hits sit tens of bp away.
3. Keep the unit-richest group.

Internal structure is preserved. A TH01 read yields
`AATGAATGAATGAATGAATGAATGATGAATGAATGAATG`, with the `ATG` spacer inside the
core. Only flanks are trimmed.

Markers whose canonical motif is on the minus strand (vWA, D5S818, CSF1PO, …)
are reverse-complemented first. A read with no detectable motif run falls back
to window length.

### 2.2 Length, not repeat-unit count

The bin key is core length in bases. Unit count is not equivalent:

| Allele | Structure | Units | Core length |
|---|---|---|---|
| TH01 9 | `[AATG]9` | 9 | 36 bp |
| TH01 9.3 | `[AATG]6 ATG [AATG]3` | 9 | 39 bp |

Both have nine units, so a unit count merges them. Core length separates them.
The same applies to any scheme that rounds lengths together.

### 2.3 Split by sequence identity

Within each bin, reads are grouped by pairwise edit-distance identity (edlib,
`--identity`, default 0.97), seed-and-grow.

This separates iso-alleles: two alleles of equal core length but different
internal structure share a bin and differ only by sequence.

### 2.4 Output

TH01 in the demo sample gives two clusters of 10 and 7 reads, genotyped
`9.3(10) / 7(7)`, plus eight singletons.

Two ONT reads of one allele diverge by 2 to 4%, near the identity threshold, so
some reads of a real allele form singleton clusters. These fall below
`--analytical-thresh` and are classified as noise. They appear in `--trace` as
candidates and are not counted toward any allele's depth.

## 3. Consensus: clusters to sequences

**Module:** `frontstr/evidence/consensus.py`

Each cluster's reads are collapsed into one consensus by partial-order
alignment (SPOA), correcting individual read errors by per-column vote.

The ISFG string, the iso-allele match, the allele number and the VCF `ALT` are
all computed from this sequence.

Alignment is global, not local. Clusters are length-binned, so their members
are equal-length full-window sequences; local alignment could trim flanks and
change the called length, which is the allele number.

### A POA backend is required

Without one the consensus falls back to the most frequent exact sequence, which
is a single raw read including its errors. On a synthetic 252 bp locus at ONT
R10 error rates:

| Backend | 3 to 4 reads | 10 to 16 reads |
|---|---|---|
| POA | 0 edits from truth | 0 edits |
| Most-common-sequence | 4 to 6 edits | 2 to 3 edits |

On five real samples the fallback produced four false microvariants in 202
called alleles:

| Sample | Marker | Fallback | POA |
|---|---|---|---|
| HG00097 | TH01 | 6.3 | 7 |
| HG00154 | TH01 | 6.1 | 6 |
| HG00097 | D13S317 | 14.1 | 14 |
| HG00263 | D18S51 | 11.3 | 12 |

TH01 6.3 and D18S51 11.3 are not alleles in the forensic literature. At
HG00263 D18S51 the unpolished consensus carried a spurious `AGAA` in the left
flank, a read error that spelled the motif and corrupted the bracket structure.

`pyspoa` is the default backend. `pyabpoa` hardcodes AVX2 x86 intrinsics and
does not build on macOS arm64. Each cluster records which method produced its
consensus, and the fallback raises `CONSENSUS_FALLBACK` on every affected
marker.

---

## 4. ISFG nomenclature and allele number

**Module:** `frontstr/interp/naming.py`, `isfg.py`, `allele_numeric.py`

Two things are derived from each consensus.

**The ISFG bracket string.** A greedy left-to-right scan taking the longest
motif run at each position, for example `[AATG]6 ATG [AATG]3`. Non-motif bases
appear in lowercase. Minus-strand markers are reverse-complemented first, so
the output is in canonical orientation.

**The allele number**, from STRNaming where available. The ISFG DNA Commission
(Gettings et al. 2024) names STRNaming as the program that produces the
designation rather than a rule set to reimplement. The number is what STRNaming
reports for the marker's standard reporting range: `CE29_TCTA[4]TCTG[6]…`
gives 29.

Two implementation details:

- The reporting range is located in the consensus by aligning the reference
  flanks just outside it, not by reference coordinates. STRNaming's ranges hug
  the array (TPOX's begins on its first base), and coordinate slicing loses the
  extra units of a long allele: it names both HG00113 TPOX alleles CE8 where
  the answer is 9 and 11.
- Reference sequence comes from a committed GRCh38 slice cache
  (`frontstr/interp/data/strnaming_ranges.tsv`, about 13 kB, built by
  `frontstr/panel/seed_strnaming.py`). Naming needs no network call and no
  `--reference`, and is reproducible byte for byte.

For markers STRNaming defines no range for (DYS393, AMEL), and for any
consensus whose range cannot be located, the number comes from arithmetic:

- **`period > 0`**, simple markers. From length:
  `divmod(length − corr_value, period)`. The integer part is full repeats, the
  remainder is the microvariant decimal. TH01 at 39 bp gives `divmod(39, 4)` =
  (9, 3), so 9.3. `corr_value` is the non-repeat content of the window,
  calibrated per marker against GRCh38.
- **`period = −1`**, compound markers (vWA, FGA, D21S11, …). By counting repeat
  units in the bracket string, minus a calibrated correction.

`Allele.number` holds the chosen number, `Allele.number_method` how it was
derived, and `Allele.number_label` the string every view renders.

---

## 5. Stutter model and classification

**Modules:** `frontstr/interp/stutter.py`, `classify.py`,
`frontstr/panel/stutter_calib.py`

Clusters above 20% of locus coverage become parents. For each parent the model
generates the sequences its stutter products would have, the longest and
second-longest uninterrupted motif runs shortened by 1 or 2 units or lengthened
by 1, with an expected read count for each.

Every cluster is classified against those expectations and against the coverage
thresholds:

| Status | Meaning |
|---|---|
| `allele` | A real allele |
| `stutter` | Matches an expected stutter product at expected depth |
| `artefact` | Above `--analytical-thresh`, below `--calling-thresh` |
| `noise` | Below `--analytical-thresh` |
| `inexact_allele` | Real, but reconstructed rather than observed |

### The model

Stutter rate is log-linear in LUS, the longest uninterrupted stretch of the
motif in the parent allele:

```
rate(−1)   = exp(−12.1125 + 0.7159 × clamp(LUS, 10, 14))
rate(step) = rate(−1) × {−1: 1.0, −2: 0.242, +1: 0.726}
```

LUS dominates. Measured over 76 observations at 52 loci on five ONT R10
samples:

| Parent LUS | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|
| Rate at −1 | 0.010 | 0.012 | 0.035 | 0.060 | 0.122 |

That is more than a tenfold range across the LUS values where most CODIS
alleles sit, so a constant rate cannot represent it: one value over-suppresses
candidates at short-LUS loci and under-suppresses at long-LUS loci at the same
time. The log-linear fit gives R² 0.965. A linear fit reaches R² 0.29 and
crosses zero near LUS 9, inside the range where stutter is still observed.

Forward stutter runs at 0.73 of the reverse rate rather than the 0.5 typical of
PCR-based protocols. On PCR-free ONT this signal is largely sequencing indel
error inside the array, which is more symmetric than polymerase slippage.

> **The model is fitted on PCR-free WGS and contains no PCR slippage
> component.** Amplicon casework stutters more than it predicts and needs a
> re-fit. `StutterModel.protocol` records the regime so the mismatch is
> visible. Re-fit with
> `frontstr calibrate-stutter --bam … --panel … --protocol …`.

Details of the calibration are in
[`stutter_calibration.md`](stutter_calibration.md).

---

## 6. Haplotype-aware phantom suppression

**Module:** `frontstr/interp/haplotype.py`

A diploid locus carries one allele per haplotype. Two candidates assigned to
the same haplotype and differing by less than one repeat unit cannot both be
real.

The weaker is demoted to `HP_PHANTOM`, its reads are recorded on the owner's
`n_reads_absorbed`, and `HP_PHANTOM_COLLAPSED` is raised.

Suppression fires only when all of these hold:

- both clusters carry at least 3 phased reads at 80% or more haplotype purity,
  so the rule is a no-op on unphased BAMs
- they share that haplotype
- they differ by 2 bp or less, well under one repeat unit
- the marker does not set `allow_triallelic`, since a duplication genuinely
  puts two alleles on one haplotype

A read-count floor cannot substitute for it. GM19038 D12S391:

```
289 bp · 10 reads · HP1        real allele
288 bp ·  6 reads · HP1        phantom of the above
301 bp ·  5 reads · HP2        the other real allele
302 bp ·  4 reads · HP2        phantom of that one
```

The two strongest clusters are both phantoms of one haplotype, and the genuine
second allele has 5 reads. Raising a read floor deletes it and reports a
homozygote. The haplotype tag resolves the locus to `289 / 301`.

With repeat-core binning in place the rule fires at 1 locus in about 110.

### Haplotype labels are local to a phase block

`HP` is meaningful only inside its own `PS`. HP1 in one block and HP1 in the
next are unrelated: the phasing tool restarted numbering, it did not switch
chromosome. Every rule keys on `(phase_set, hp)`, and a cluster whose tagged
reads span more than one block is treated as unphased however pure its HP
labels look.

Across the five ONT slices, 3 of 125 loci have spanning phased reads from more
than one block. HG00097 D13S317 is one: a 14-read cluster reporting 100% HP2,
drawn 4 reads from one block and 3 from another.

Those loci raise `PHASE_BLOCK_SPLIT` (WARN) and `--trace` marks the cluster
`[2 phase blocks, HP not comparable]`.

A BAM with `HP` but no `PS` keeps the unscoped behaviour, so haplotype
reasoning is not disabled for phasing tools that omit the tag.

### The same invariant in reverse

If two candidates on the same haplotype cannot both be real, two on opposite
haplotypes cannot be the same allele, whatever their read ratio.
`on_opposite_haplotypes()` exposes that, and stage 8 uses it to stop the
peak-height ratio collapsing a genuine heterozygote. Same evidence and same
gates, opposite direction.

---

## 7. Iso-allele catalog (optional)

**Module:** `frontstr/interp/catalog.py`

With `--catalog`, each called allele's repeat core is compared by edit distance
against curated sequences. An exact match adopts the published ISFG string and
the iso-allele suffix. Within 2 edits the match is recorded as approximate.
Beyond that the live-computed values stand.

This is what distinguishes D3S1358 allele `14` from `14b`: the same repeat
count with a distinguishable sequence.

The shipped catalog (`examples/catalogs/demo_seed.json`) is a four-entry
demonstration seed, not an STRSeq import.

---

## 8. Genotype call

**Module:** `frontstr/interp/triallelic.py`

Surviving candidates, sorted by read count, become a genotype:

| Candidates | Outcome |
|---|---|
| 0 | `no_data` |
| 1, or 2nd below `--min-phr` of the 1st and not on the opposite haplotype | `homozygous` |
| 2 | `heterozygous` |
| 3 or more, marker allows triallelic | `triallelic_type_I` / `type_II` / `review`, by balance |
| 3 or more, marker does not | `heterozygous`, or `mixture_suspected` if the 3rd is substantial |

A 3rd candidate must clear both `--calling-thresh` and the absolute floor
`--min-reads-third`.

### Haplotype evidence overrides the peak-height ratio

`--min-phr` is inherited from capillary electrophoresis, where peak height is
the only evidence. On phased long reads it is not. HG00113 D2S1338:

```
20 · 17 reads · 100% HP1
17 ·  5 reads · 100% HP2     PHR 0.29, below the 0.4 floor
```

On read ratio alone that is homozygous `20`. Illumina, and two other callers,
report `17/20`. At ONT depths a 5-against-17 split is ordinary sampling, and
the phasing assigns the two candidates to opposite haplotypes.

The ratio therefore yields to `on_opposite_haplotypes()` (stage 6), under the
same gates: at least 3 tagged reads at 80% purity on both sides, so it is a
no-op on unphased BAMs. The rescued allele carries `hp_rescued` and the marker
raises `HP_RESCUED_HET`.

Across the five test samples the rescue fires at one locus.

**Amelogenin takes a separate path** (`interp/amel.py`), since it is not a
tandem repeat. It counts reads at the AMELX and AMELY regions and reports X, Y
or both.

---

## 9. QC flags

**Module:** `frontstr/interp/qc.py`

Run-level conditions, applied after every marker is called:

| Flag | Fires when | Severity |
|---|---|---|
| `DROPOUT` | No allele called at the locus | warn |
| `LOW_COVERAGE` | Called below `--low-coverage-reads` | warn |
| `STRAND_BIAS` | An allele's strand ratio is skewed beyond chance | warn |
| `ALLELE_IMBALANCE` | A heterozygote's balance exceeds `--balanced-ab-max` | warn |
| `INEXACT_ALLELE` | A called allele is a reconstruction | info |
| `CE_NOMENCLATURE_OFFSET` | The number is knowingly not the kit designation | warn |

Raised earlier in the pipeline: `MIXTURE_SUSPECTED`, `TRIALLELIC`, `ISOALLELE`,
`CONSENSUS_FALLBACK`, `HP_PHANTOM_COLLAPSED`, `HP_RESCUED_HET`,
`PHASE_BLOCK_SPLIT`, `NON_DEFAULT_THRESHOLD`.

### The coverage floor

`--low-coverage-reads` is measured against the reads supporting the genotype,
not against every read spanning the window. Reads that clustered into neither
allele are not draws from the pair the derivation models. Measuring the
spanning total instead let HG00263 D18S51 pass unflagged: called on 11 reads of
33 spanning, missing the second allele Illumina reports.

The risk being modelled is a true heterozygote whose minor allele falls below
`--calling-thresh` and is not called, reporting a homozygote. Taking the most
unbalanced heterozygote the caller accepts (`--min-phr` 0.4, so the minor
allele is 28.6% of the pair):

| Floor | Dropout risk | Loci flagged |
|---|---|---|
| 17 | 9.7% | 12% |
| **20** | **5.7%** | **25%** |
| 25 | 4.6% | 43% |

20 is the knee. Moving 17 to 20 halves the risk for 13 points of flag rate;
moving 20 to 25 buys one point of risk for eighteen.

A quarter of loci flagged at about 30x ONT is the expected rate. Median called
coverage in the reference slices is 26 and the lower quartile is 20, so the
floor sits at Q1.

### Strand bias

An exact two-sided binomial test. Alleles below 10 reads are not tested,
because the test cannot reach `strand_bias_p` even for a perfect 0/n split. On
151 real called alleles it fires once, about the rate chance predicts at
p < 0.01.

### Nomenclature offset

`CE_NOMENCLATURE_OFFSET` is curated in the panel through
`System.kit_nomenclature_note`, set for vWA and D21S11. Which markers diverge
from a kit convention is a property of the kit a laboratory compares against.

The flag is suppressed when every called allele was named by STRNaming, since
the number is then the kit designation. The notes remain in the panel because
they describe the fallback path.

---

## 10. Elsewhere

| Topic | Location |
|---|---|
| Installing and running, commands, flags, output formats | [`../README.md`](../README.md) |
| The ten parameters, defaults and provenance | [`../README.md`](../README.md), and `frontstr/params.py` |
| How each default was derived | the docstring beside it, mostly `frontstr/interp/qc.py` |
| Stutter calibration | [`stutter_calibration.md`](stutter_calibration.md) |
| Module layout | [`architecture.md`](architecture.md) |
| The audit trail | `frontstr/audit.py` |
| Demo sample and its expected profile | [`../demodata/README.md`](../demodata/README.md) |
| Known limitations, planned work | [`../README.md`](../README.md), [`../ROADMAP.md`](../ROADMAP.md) |
