# The ten parameters

Every knob a FRONTStr run can turn. For each one: what it does, the exact line
of code that applies it, where its default came from, and what changing it
costs in both directions.

The table lives in [`frontstr/params.py`](../frontstr/params.py) and is printed
at the start of every run. That is deliberate: the values that decide a call are
usually the ones nobody typed, so a run listing only its command line reads as
if it were barely configured.

- [How to read this](#how-to-read-this)
- [Reads kept](#1-reads-kept) — `min_mapq`, `flank_anchor`
- [Grouping reads](#2-grouping-reads) — `identity_threshold`, `len_tolerance_bp`
- [What counts as an allele](#3-what-counts-as-an-allele) — `analytical_thresh`, `calling_thresh`
- [The genotype call](#4-the-genotype-call) — `min_phr_for_het`, `min_reads_third`
- [QC policy](#5-qc-policy) — `low_coverage_reads`, `balanced_ab_max`
- [Not in the table](#not-in-the-table)
- [References](#references)

---

## How to read this

### Provenance is part of the value

Each parameter carries where its default came from. This is not decoration: it
decides what happens when you override it.

| Provenance | Meaning | Overriding it |
|---|---|---|
| `derived` | Computed from measured data | **Marks the run.** Raises `NON_DEFAULT_THRESHOLD` on every marker |
| `chosen` | Defensible but picked. Nothing measured says it must be this | Ordinary tuning |
| `convention` | Inherited from forensic practice, not from this data | Ordinary tuning |

Three of the ten are `derived`: `len_tolerance_bp`, `min_reads_third`,
`low_coverage_reads`. Overriding any of them raises the flag at **marker**
level, not run level, because that is where every consumer already looks: the
audit census, the XLSX QC sheet, the HTML row tint. A run-level note is the one
thing a reviewer scanning per-locus rows never sees.

The point is not to stop you testing a threshold. `--min-reads-third 2` is a
legitimate experiment. It is that six months later the profile still says which
numbers were changed, and which of them had measured backing.

### Where a value ends up

One override reaches four places:

1. The **parameter block** printed at the start of the run, with `CHANGED`
   beside it and the default it replaced.
2. The **`AuditRecord`** embedded in the canonical JSON, via
   `RunParameters.as_audit_rows()`, so the record travels with the result.
3. The **run header** in `--trace`, which lists overrides and calls out the
   `derived` ones separately.
4. `NON_DEFAULT_THRESHOLD` on every marker, if the override was `derived`.

`tests/test_params.py` asserts the table's defaults equal the ones the pipeline
actually uses. A table that drifts is worse than none, because the report would
confidently print a default nothing is using.

### A note on what is *not* here

The **stutter model** is not a parameter. It is a calibrated object with its own
version, provenance and fit statistics, and it is named in the trace header
rather than in this table. See [stutter_calibration.md](stutter_calibration.md).

---

## 1. Reads kept

Applied in [`frontstr/evidence/pileup.py`](../frontstr/evidence/pileup.py), in
`_read_to_observation`, in this order. The first failing test rejects the read
and its reason is counted, which is what the `Rejected` block of `--trace`
enumerates.

### `min_mapq` — 20 — *chosen*

**What it does.** Drops reads whose aligner-assigned mapping quality is below
the threshold.

```python
if read.mapping_quality < min_mapq:
    return None, RejectReason.LOW_MAPQ
```

**Where the value comes from.** Convention, not measurement. MAPQ 20 is the
long-standing "roughly 99% confident this read is placed correctly" line
(MAPQ = −10 log₁₀ P(wrong placement), so 20 → P = 0.01), and it is the default
in most short-variant pipelines. Nothing in FRONTStr's own data was used to pick
it.

**Why it matters more here than elsewhere.** This is also what keeps X/Y
paralogue mismappings out of the sex markers. AMELX and AMELY are homologous;
reads from one can align to the other with low MAPQ. Dropping the ambiguous ones
is what makes the X/Y read counts mean anything.

**Raising it** costs coverage in repetitive flanks, which is where CODIS loci
sit. Depth is already the binding constraint on ONT panels, so this is not free.

**Lowering it** admits reads whose placement is uncertain. In a repeat region an
ambiguously placed read is more likely to carry a *different* locus's repeat
structure, which does not look like noise once it clusters. It looks like an
allele.

### `flank_anchor` — 20 — *chosen*

**What it does.** Requires the read to extend at least this many bases past the
marker window on **both** sides before it can become an observation.

```python
if ref_start > start - flank_anchor:
    return None, RejectReason.LEFT_FLANK_SHORT
if ref_end < end + flank_anchor:
    return None, RejectReason.RIGHT_FLANK_SHORT
```

**Why it exists at all.** A read that stops inside the repeat has a *truncated*
repeat tract, and the truncation is indistinguishable from a shorter allele.
Admitting one is not adding noise, it is adding a fake allele of a plausible
length. Requiring clean flank on both sides means every surviving read spans the
whole array, so its length is the allele's length.

**Where 20 comes from.** Chosen. It is enough sequence to anchor the alignment
outside the repeat, and short enough not to lose reads at ONT read lengths where
the panel windows are already widened to the repeat boundary ±100 bp.

**Interaction with the panel.** The panel windows are wider than the repeat by
design, so `flank_anchor` is measured from the *window* edge, not the repeat
edge. The effective anchor outside the repeat is therefore ~120 bp.

**Raising it** rejects more reads and is the fastest way to make a locus
`no_data`. **Lowering it** re-admits partially spanning reads, whose truncated
tracts bin as shorter alleles.

> When a locus reports no usable reads, `--trace` names these two parameters
> explicitly, because between them they account for nearly every rejection.

---

## 2. Grouping reads

Applied in [`frontstr/evidence/cluster.py`](../frontstr/evidence/cluster.py).
Reads are first binned by repeat-core length, then each bin is split by pairwise
sequence identity.

### `identity_threshold` — 0.97 — *chosen* ⚠️ known to be miscalibrated

**What it does.** Pairwise identity required for a read to join an existing
cluster's seed.

```python
# _identity(a, b) = 1 - editDistance / max(len(a), len(b)), Levenshtein via edlib
if _identity(seed.sequence, m.sequence) >= identity_threshold:
    cluster_members.append(m)
```

**The problem, stated plainly in the code.** The value was set as if comparing a
read to a *consensus*. The comparison is actually **raw read to raw read**, and
two ONT reads of the same allele diverge by 2–4%. At 0.97 the threshold sits
right on that boundary, so genuine same-allele reads are split into separate
clusters.

**What that does downstream.** Splitting inflates the cluster count. Most of the
fragments fall under `analytical_thresh` and are discarded as noise, so the
genotype usually survives, but the fragments show up as candidates in the trace
and they eat into `called_reads`. This is visible in the example data: HG00113
D3S1358 produces 5 clusters from 4 length bins, four of them single reads.

**Why it has not been changed.** Changing it moves every locus at once, and the
right value has to be measured rather than guessed. It is listed here as a known
open issue, not as a tuned value. It is the parameter most worth experimenting
with.

**Raising it** splits more. **Lowering it** merges more, and past some point it
will merge a genuine microvariant into its neighbour — the failure mode is a
lost allele, not a spurious one.

### `len_tolerance_bp` — 0 — *derived* — **do not change**

**What it does.** Merges adjacent length bins whose keys differ by at most this
many bases, before identity clustering.

```python
while j < len(keys) and keys[j] - anchor <= len_tolerance_bp:
    group.extend(bins[keys[j]])
```

**Why it is zero.** It was introduced to absorb indel errors in the flanks. That
problem no longer exists: reads are binned by **repeat-core length**, not by
window length, so flank indels never reach the binning key.

**Why it must stay zero.** With the flank problem gone, any tolerance merges
genuine microvariants. The canonical case is TH01 **9** and **9.3**, which are
3 bp apart. A tolerance of 3 collapses them into one allele. That is a false
homozygote at one of the most-used loci in forensic genetics.

**This is the parameter where a well-meaning override does the most damage**,
which is why it is `derived`: changing it marks every marker in the run.

---

## 3. What counts as an allele

Applied in [`frontstr/interp/classify.py`](../frontstr/interp/classify.py). Both
are fractions of the locus's **total spanning reads**, not of the called
coverage. This is the one place the spanning total is the right denominator:
it is what decides where the noise line falls, and changing it would move what
counts as noise.

```python
es = expected_stutter.get(allele.consensus, 0.0)
if es > 0 and allele.n_reads_total <= es:
    return AlleleStatus.STUTTER

frac = allele.n_reads_total / total_reads
if frac < analytical_thresh:
    return AlleleStatus.NOISE
if frac < calling_thresh:
    return AlleleStatus.ARTEFACT
return AlleleStatus.ALLELE
```

Note the order: the **stutter model runs first**. A candidate at a stutter
position with a count at or below the model's expectation is stutter regardless
of the two thresholds. Only what survives that is measured against them.

### `analytical_thresh` — 0.02 (2%) — *chosen*

**What it does.** Below this fraction of locus coverage, a cluster is `NOISE`
and disappears from the report entirely.

**Where it comes from.** The analytical threshold is a standard concept in
forensic DNA interpretation: the level below which a signal is not
distinguishable from baseline. The 2% figure is FRONTStr's own choice; it is not
transferred from a CE validation, where the analytical threshold is in RFU and
means something physically different.

**Raising it** hides low-level candidates, including real minor alleles.
**Lowering it** fills the trace with basecaller noise. Because it only decides
*visibility*, not calling — that is `calling_thresh` — a lower value is the
safer direction when investigating a locus.

### `calling_thresh` — 0.10 (10%) — *chosen*

**What it does.** Above `analytical_thresh` but below this, a cluster is
`ARTEFACT`: real enough to show, not real enough to genotype.

**Where it comes from.** Chosen, and the code says so in as many words. It is
not data-derived. 10% is a common stochastic-threshold order of magnitude in
forensic practice, but no measurement in this project fixes it there.

**It is load-bearing well beyond classification.** It appears in the
`low_coverage_reads` derivation below, and in `call_profile` when deciding
whether a third candidate is real. Moving it moves both.

**Raising it** collapses unbalanced heterozygotes into homozygotes — a false
exclusion, the most serious error class in forensic identification.
**Lowering it** admits artefacts as alleles, tending toward false triallelic
calls and false mixture flags.

---

## 4. The genotype call

Applied in [`frontstr/interp/triallelic.py`](../frontstr/interp/triallelic.py),
`call_profile`.

### `min_phr_for_het` — 0.4 — *convention*

**What it does.** The second candidate must reach this fraction of the first, by
read count, or the locus is collapsed to homozygous.

```python
if phr_12 < min_phr_for_het:
    # collapse to homozygous
```

**Where it comes from.** Capillary electrophoresis, where **peak height is all
there is**. The 60/40 heterozygote balance expectation is long-standing forensic
practice, and 0.4 is that convention transcribed to read counts.

**Why the provenance matters here more than anywhere else.** On phased long
reads, peak height is *not* all there is. A diploid locus carries one allele per
haplotype, so two candidates confidently on opposite haplotypes are two alleles
whatever their ratio. At ONT depths a 5-vs-17 split (PHR 0.29) is ordinary
sampling, and collapsing it is a false exclusion waiting to happen.

So this floor is **overridden by haplotype evidence**:
`frontstr.interp.haplotype.on_opposite_haplotypes` rescues the second allele and
marks it `hp_rescued`, raising `HP_RESCUED_HET`. It is the same invariant that
`suppress_hp_phantoms` uses to *delete* alleles, applied in the opposite
direction. Both are no-ops on an unphased BAM, where the ratio is again the only
evidence there is.

On the reference slices the rescue fires at exactly one locus. It is meant to
stay that narrow.

**Raising it** collapses more heterozygotes. **Lowering it** promotes stutter
and artefacts into second alleles.

### `min_reads_third` — 5 — *derived*

**What it does.** An absolute read floor, **in addition** to the fractional
`calling_thresh`, that a third candidate must clear before it can promote a
locus to triallelic or raise `MIXTURE_SUSPECTED`.

**Where it comes from.** The known-bug #6 work on ONT basecaller phantoms.
Measured behaviour: a phantom of **2–4 reads sits at 5–10% allele fraction** at
the 20–50× coverage these panels run at. That is *above* `calling_thresh` at the
low end of the coverage range, so a fraction alone cannot exclude it. An
absolute floor can.

**Why a fraction was not enough.** The two thresholds fail in opposite regimes.
At 20× a 3-read phantom is 15% and clears any sane fraction; at 50× a genuine
minor contributor at 8% is 4 reads and is excluded by any sane absolute floor.
FRONTStr requires **both**, which is conservative in the direction that matters:
a false mixture flag on a single-source sample is a serious error.

**Lowering it re-admits the phantoms**, which is why it is `derived`. A
per-marker override exists in the panel (`System.min_reads_third`) for loci with
genuine triallelic propensity.

---

## 5. QC policy

Applied in [`frontstr/interp/qc.py`](../frontstr/interp/qc.py). These two do not
change any genotype. They decide what gets flagged.

### `low_coverage_reads` — 20 — *derived*

**What it does.** A called locus supported by fewer reads than this raises
`LOW_COVERAGE`.

```python
if result.called_reads < thr.low_coverage_reads:
```

**Measured against `called_reads`, not spanning reads.** Reads that clustered
into neither allele are not draws from the pair the derivation models. This is
not theoretical: using the spanning total hid a real dropout. HG00263 D18S51 is
called on 11 reads out of 33 spanning, misses the second allele Illumina sees,
and raised no flag.

**The derivation, re-run against ONT data (2026-07).** The risk being modelled
is a true heterozygote whose minor allele falls under `calling_thresh` and is
therefore not called, so a homozygote is reported: a false exclusion.

1. Take the most unbalanced heterozygote the caller still accepts. With
   `min_phr_for_het` = 0.4, the minor allele is 0.4 / 1.4 = **28.6%** of the
   pair.
2. Ask how often, by binomial sampling, it lands below `calling_thresh` × the
   spanning total.
3. Convert spanning to called coverage: called is a median **0.795** of spanning
   across 117 called loci in the reference slices.

The curve is a **sawtooth**, because the floor is an integer and the risk jumps
whenever it crosses one. Read as "from this coverage upward the risk never again
exceeds":

| Floor | Residual risk | Loci flagged |
|---|---|---|
| 12 | 12.2% | 3% |
| 15 | 12.2% | 8% |
| 17 | 9.7% | 12% |
| **20** | **5.7%** | **25%** ← the knee |
| 22 | 5.7% | 31% |
| 25 | 4.6% | 43% |

20 is where the curve turns. 17 → 20 halves the risk for 13 points of flag rate;
20 → 25 buys one point of risk for eighteen.

**A correction this derivation forced.** The docstring previously claimed ~1.1%
risk at N = 20. That figure came from measuring against the spanning total and
is wrong under the corrected model. It is **~5.7%**.

**A quarter of loci flagged at ~30× ONT is the expected rate, not a symptom.**
Median called coverage in the reference slices is 26 and the lower quartile is
20, so the floor sits at Q1 by construction.

### `balanced_ab_max` — 0.65 — *chosen*

**What it does.** The largest allele balance a heterozygote may have and still
count as balanced. Above it, `ALLELE_IMBALANCE`.

**The scale matters, and it is one-sided.** Allele balance here is the
**strongest called allele over the sum of the called pair**, so it runs from
0.50 (perfectly even) to 1.0 (everything on one allele). It is not symmetric
around 0.5, and a band written as if it were would be wrong.

**Where 0.65 sits.** `min_phr_for_het` = 0.4 corresponds to **0.714** on this
scale, the point below which no heterozygote is called on read counts at all. So
0.65 sits *inside* the callable range: it warns about heterozygotes that were
called but are close to the edge, rather than rejecting anything.

**Raising it toward 0.714** makes the flag fire only just before the call would
have been collapsed. **Lowering it toward 0.5** floods the report; at ONT depths
ordinary sampling produces balances in the high 0.5s routinely.

---

## Not in the table

Two QC thresholds live on `QcThresholds` and are not exposed on the CLI:

- **`strand_bias_p`** = 0.01. Two-sided binomial p-value below which a strand
  ratio is called biased. Kept strict deliberately: at ONT panel coverages a 5%
  cutoff fires on ordinary sampling noise often enough to train reviewers to
  ignore the flag.
- **`strand_bias_min_reads`** = 10. Alleles below this are not strand-tested at
  all. Under 10 reads the exact binomial cannot reach `strand_bias_p` even for a
  perfect 0/n split, so testing would only ever produce false reassurance.

Panel-level values (`corr_value`, `period`, `reference_ce`, `strand`,
`allow_triallelic`, per-marker `min_reads_third` and `ont_len_tolerance`) are
properties of a marker, not of a run. They live in the panel YAML.

---

## References

Only sources this project actually used are listed. Where a default was picked
rather than derived, this document says so instead of attaching a citation to
it.

- **Gettings et al. 2024**, ISFG DNA Commission recommendations on sequence
  data. Recommendation 2 names STRNaming as the program that produces bracketed
  repeat formatting; this is why FRONTStr takes both the allele number and the
  bracket string from STRNaming rather than computing its own. See
  [`frontstr/interp/naming.py`](../frontstr/interp/naming.py) and
  [how_it_works.md §4](how_it_works.md).
- **[docs/stutter_calibration.md](stutter_calibration.md)** — the measured ONT
  stutter model, its fit, and how to re-fit it. The model is what runs *before*
  `analytical_thresh` and `calling_thresh`.
- **[docs/how_it_works.md](how_it_works.md)** — the pipeline these parameters
  act on, step by step. §2 clustering, §5 stutter, §8 the genotype call,
  §9 QC flags.
- **Known-bug #6 (ONT basecaller phantoms)** — the measurement behind
  `min_reads_third`. Internal to this project; the numbers are in
  [`frontstr/panel/models.py`](../frontstr/panel/models.py).
- **MAPQ** as −10 log₁₀ P(wrong placement) is the SAM specification's
  definition; the choice of 20 is conventional, not derived here.

**Not cited, on purpose.** `analytical_thresh`, `calling_thresh` and
`min_phr_for_het` have recognisable counterparts in CE validation practice, but
the numbers here were not transferred from a published CE validation and
attaching one would overstate their backing. `min_phr_for_het` is the closest to
inherited, and its provenance is recorded as `convention` for exactly that
reason.

---

## Seeing the values for a run

Every run prints them. To see them without running anything:

```bash
frontstr interpret --bam SAMPLE.bam -p examples/panels/codis_20_grch38.yaml --show-params
```

To watch one act on a locus, including which candidates each threshold
discarded and why:

```bash
frontstr interpret --bam SAMPLE.bam -p examples/panels/codis_20_grch38.yaml --trace
```

To override one and see the run mark itself:

```bash
frontstr interpret --bam SAMPLE.bam -p examples/panels/codis_20_grch38.yaml --min-reads-third 2 --trace
```

The header will list the override, call out that a `derived` default was
changed, and every marker will carry `NON_DEFAULT_THRESHOLD`.
