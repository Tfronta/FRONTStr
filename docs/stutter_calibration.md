# Stutter calibration — ONT R10, first pass

**Status:** first pass, deliberately provisional. Fitted on 5 samples.
**Model shipped as:** `frontstr.panel.stutter_calib.DEFAULT_STUTTER_MODEL`, version `2026.07-ont-r10-wgs`.

> **The one thing to read if you read nothing else:** this model is fitted on
> **PCR-free WGS**. There is no PCR slippage component in it, because there was
> no PCR. Casework on an amplicon panel will have more stutter than this model
> predicts, and the model must be re-fitted on amplicon data before it is used
> for casework. `StutterModel.protocol` records this so it cannot be lost.

---

## 1. What was wrong before

FRONTStr shipped with the stutter constants inherited from toaSTR, which are CE
/ Illumina figures:

| Parameter | Old value | Meaning |
|---|---|---|
| `default_lus` | 0.10 | per-step rate inside the longest uninterrupted run |
| `default_slus` | 0.05 | same for the second-longest run |
| `plus_factor` | 0.5 | forward stutter relative to reverse |
| step decay | `rate ** k` | geometric: −2 is the −1 rate squared |

Three things are wrong with this on ONT, and the third is the serious one.

## 2. Measurement

Method, implemented in `frontstr/panel/stutter_calib.py`:

- Only loci where a stutter position **cannot** be confused with a real allele:
  homozygotes, or heterozygotes at least 3 repeat units apart. On the CODIS 20
  panel this discards roughly half the loci (52 usable, 60 discarded).
- Reads grouped by **repeat-core length**, and all reads at a given offset
  summed. A stutter peak is a *position*, not a cluster — a peak split across
  two clusters must not halve the measured ratio.
- Parents with fewer than 8 reads skipped.
- **Zero-stutter positions recorded as zeros.** This is the one that bites:
  measuring only the positions where a stutter cluster exists conditions the
  estimate on stutter being present. Doing that gives 0.098 for the −1 step;
  including the 46 zero observations gives the correct **0.044**. An early
  version of this analysis made exactly that mistake and nearly shipped a rate
  more than twice too high.

Data: 5 ONT R10.4.1/LSK114 Dorado 1000G WGS slices (HG00097, HG00113, HG00154,
HG00263, GM19038), CODIS 20 + sex panel. 76 parent observations over 52 loci.

## 3. Results

### 3.1 Pooled rates

| Step | Old model predicts | **Measured** | Zero observations |
|---|---|---|---|
| −1 | 0.100 | **0.044** | 46 / 76 |
| −2 | 0.010 (= 0.10²) | **0.011** | 62 / 76 |
| +1 | 0.050 (= 0.10 × 0.5) | **0.032** | 44 / 76 |

Two findings:

- The flat −1 rate is **2.3× too high** on average.
- Forward stutter runs at **0.73** of the reverse rate, not 0.5. On PCR-free
  ONT this signal is largely sequencing error inside the repeat array, and
  indel error is far more symmetric than polymerase slippage. The −2 step is
  0.24 of the −1 rate rather than the 0.10 a geometric decay implies — again
  because this is not a multi-cycle slippage process.

### 3.2 The dominant effect: LUS

| LUS | n | Pooled −1 ratio | Fitted |
|---:|---:|---:|---:|
| 10 | 7 | 0.0100 | 0.0070 |
| 11 | 20 | 0.0121 | 0.0145 |
| 12 | 14 | 0.0346 | 0.0296 |
| 13 | 10 | 0.0604 | 0.0605 |
| 14 | 10 | 0.1222 | 0.1237 |

*(bins with n < 7 — LUS 4, 6, 7, 8, 9, 15, 16, 17 — are excluded from the fit;
see §5)*

The rate spans **more than 10×** across the range where most CODIS alleles sit.
No single constant can represent that. A flat 0.10 is simultaneously ~10× too
high at LUS 10 and slightly too low at LUS 14 — so the old model was
over-suppressing candidates at short-LUS loci and under-suppressing at long-LUS
ones, which is the worst of both.

### 3.3 Model form

```
rate(-1)   = exp(-12.1125 + 0.7159 × clamp(LUS, 10, 14))     R² = 0.965
rate(step) = rate(-1) × {-1: 1.0, -2: 0.242, +1: 0.726}
```

A plain linear fit was tried first and **rejected**: the measured rates
accelerate with LUS, so a straight line fits poorly (R² = 0.29 over all bins,
0.52 over hand-picked ones) and — worse — crosses zero around LUS 9, inside the
range where stutter is still observed. That would mean *no stutter model at
all* for short-LUS loci. The log-linear form is convex by construction, cannot
go negative, and treats the rate as a multiplicative per-unit hazard, which is
what a slippage process is. R² goes from 0.29 to **0.965**.

Outside LUS 10–14 the LUS is **clamped, not extrapolated**. Nothing measured
here says anything about LUS 4 or LUS 30.

### 3.4 SLUS

Not separately calibrable from this measurement: the −1 position is the same
whichever run slipped, so the LUS and SLUS contributions are not separable.

The model handles it structurally instead — the rate is a function of *the
slipping run's own length*, so a shorter secondary run automatically gets a
lower rate. The old flat `slus = lus / 2` was a crude stand-in for exactly
that. `slus_factor` remains as a tuning knob, defaulting to 1.0, and is
explicitly **not** derived from data.

## 4. Which way we chose to be wrong

Over-predicting stutter classifies a real minor allele as stutter and removes
it — a silent false homozygote the analyst never sees. Under-predicting lets a
stutter peak survive as an allele candidate, which surfaces as a review flag
and is caught downstream by `min_reads_third` and haplotype-aware suppression.

The second failure mode is far safer, so where the data is thin the model is
deliberately left conservative. Verified: HG00113's reference profile is
unchanged across all 25 markers under the new model, and `mixture_suspected`
across the 5 samples stays at 0.

## 5. Known limitations — read before trusting this

1. **PCR-free only.** Stated above; restated here because it is the limitation
   most likely to cause harm.
2. **5 samples.** 76 observations over 52 loci. Fine for establishing the shape
   of the model, thin for the parameters.
3. **Calibrated LUS range is only 10–14.** Everything outside is clamped. LUS
   8 shows a measured 0.045 (n=6, just under the inclusion threshold) which the
   clamped model under-predicts as 0.007. Widening this range is the single
   biggest reason to calibrate on more samples.
4. **No per-marker rates.** Marker-level spread is real (D18S51 0.196 vs
   D7S820 0.012 pooled) but n = 1–8 per marker is far too thin to set
   per-marker overrides. `System.stutter_overrides` exists for labs that have
   the data.
5. **Half the loci are unusable** for calibration by construction, since
   heterozygotes closer than 3 repeat units have ambiguous stutter positions.
   More samples is the only fix.
6. **The analytical/calling thresholds (0.02 / 0.10) are still not derived from
   data.** They were not part of this pass.

## 6. Next step: widen the sample set

The tooling makes this one command. Selection rule per
`feedback_ont_r10_dorado_only`: R10 + Dorado only, never R9, never guppy.

```bash
frontstr calibrate-stutter \
    --panel examples/panels/codis_20_grch38.yaml \
    --protocol wgs_pcr_free \
    --bam sample1.bam --bam sample2.bam ... \
    --out examples/stutter/ont_r10_wgs.json
```

The command prints per-LUS support and dims bins with n < 7, so a thin fit is
visible rather than silently shipped. What to watch as samples are added:

- Does the LUS range with usable support widen past 10–14?
- Does `log_slope` stay near 0.72, or drift?
- Do the step factors stay near 0.24 / 0.73?
- Does any marker separate far enough from the pooled curve to justify a
  per-marker override?
