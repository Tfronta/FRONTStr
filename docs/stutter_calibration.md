# Stutter calibration: ONT R10

**Shipped as** `frontstr.panel.stutter_calib.DEFAULT_STUTTER_MODEL`, version
`2026.07-ont-r10-wgs`. First pass, fitted on 5 samples.

> **The model is fitted on PCR-free WGS.** It contains no PCR slippage
> component, because there was no PCR. An amplicon panel stutters more than this
> model predicts and needs a re-fit before casework use.
> `StutterModel.protocol` records the regime, so a mismatch is visible in the
> audit record of any run.

---

## 1. Method

Implemented in `frontstr/panel/stutter_calib.py`.

- **Only loci where a stutter position cannot be confused with a real allele**:
  homozygotes, or heterozygotes at least 3 repeat units apart. On the CODIS 20
  panel this leaves 52 usable loci and discards 60.
- **Reads grouped by repeat-core length**, with all reads at a given offset
  summed. A stutter peak is a position, not a cluster: a peak split across two
  clusters must not halve the measured ratio.
- **Parents below 8 reads are skipped.**
- **Zero-stutter positions are recorded as zeros.** Measuring only positions
  where a stutter cluster exists conditions the estimate on stutter being
  present. On this data that gives 0.098 for the −1 step; including the 46 zero
  observations gives 0.044.

Data: 5 ONT R10.4.1 / LSK114 Dorado 1000 Genomes WGS slices (HG00097, HG00113,
HG00154, HG00263, GM19038), CODIS 20 + sex panel. 76 parent observations over
52 loci.

---

## 2. Results

### 2.1 Pooled rates

| Step | Measured rate | Zero observations |
|---|---|---|
| −1 | 0.044 | 46 / 76 |
| −2 | 0.011 | 62 / 76 |
| +1 | 0.032 | 44 / 76 |

Forward stutter runs at **0.73** of the reverse rate. On PCR-free ONT this
signal is largely sequencing indel error inside the repeat array, which is more
symmetric than polymerase slippage. The −2 step is **0.24** of the −1 rate,
rather than the 0.10 a geometric decay would give, for the same reason: this is
not a multi-cycle slippage process.

### 2.2 LUS dominates

| LUS | n | Pooled −1 ratio | Fitted |
|---:|---:|---:|---:|
| 10 | 7 | 0.0100 | 0.0070 |
| 11 | 20 | 0.0121 | 0.0145 |
| 12 | 14 | 0.0346 | 0.0296 |
| 13 | 10 | 0.0604 | 0.0605 |
| 14 | 10 | 0.1222 | 0.1237 |

Bins with n below 7 (LUS 4, 6, 7, 8, 9, 15, 16, 17) are excluded from the fit.

The rate spans more than tenfold across the range where most CODIS alleles sit,
so no single constant represents it. A flat rate is simultaneously too high at
LUS 10 and too low at LUS 14, over-suppressing candidates at short-LUS loci and
under-suppressing at long-LUS loci at the same time.

### 2.3 Model form

```
rate(-1)   = exp(-12.1125 + 0.7159 × clamp(LUS, 10, 14))     R² = 0.965
rate(step) = rate(-1) × {-1: 1.0, -2: 0.242, +1: 0.726}
```

A linear fit reaches R² 0.29 over all bins and crosses zero near LUS 9, inside
the range where stutter is still observed, which would leave short-LUS loci
with no model at all. The log-linear form is convex, cannot go negative, and
treats the rate as a multiplicative per-unit hazard.

Outside LUS 10 to 14 the LUS is **clamped, not extrapolated**. Nothing measured
here describes LUS 4 or LUS 30.

### 2.4 SLUS

The second-longest run is not separately calibrable from this measurement: the
−1 position is the same whichever run slipped, so the two contributions are not
separable.

The model handles it structurally. The rate is a function of the slipping run's
own length, so a shorter secondary run gets a lower rate automatically.
`slus_factor` remains a tuning knob, defaults to 1.0, and is not derived from
data.

---

## 3. Direction of error

Over-predicting stutter classifies a real minor allele as stutter and removes
it: a false homozygote with nothing on screen to show for it. Under-predicting
lets a stutter peak survive as a candidate, which raises a review flag and is
caught by `--min-reads-third` and haplotype-aware suppression.

The second failure mode is the safer one, so where data is thin the model is
left conservative.

Under this model HG00113's reference profile is unchanged across all 25
markers, and `mixture_suspected` across the 5 samples stays at 0.

---

## 4. Limitations

1. **PCR-free only.** Restated here because it is the limitation most likely to
   cause harm.
2. **5 samples**, 76 observations over 52 loci. Enough to establish the shape of
   the model, thin for its parameters.
3. **The calibrated LUS range is 10 to 14.** Everything outside is clamped.
   LUS 8 shows a measured 0.045 (n=6, just under the inclusion threshold) which
   the clamped model under-predicts as 0.007. Widening this range is the main
   reason to calibrate on more samples.
4. **No per-marker rates.** Marker-level spread is real (D18S51 0.196 against
   D7S820 0.012 pooled) but n = 1 to 8 per marker is too thin to set overrides.
   `System.stutter_overrides` exists for laboratories that have the data.
5. **Half the loci are unusable** for calibration by construction, since
   heterozygotes closer than 3 repeat units have ambiguous stutter positions.
6. **`--analytical-thresh` and `--calling-thresh` are not derived from data.**
   They were not part of this calibration.

---

## 5. Re-fitting

Only R10 chemistry with Dorado basecalling. R9 and guppy data will run but is
being fitted to a different error process.

```bash
frontstr calibrate-stutter \
    --panel examples/panels/codis_20_grch38.yaml \
    --protocol wgs_pcr_free \
    --bam sample1.bam --bam sample2.bam ... \
    --out examples/stutter/ont_r10_wgs.json
```

The command prints per-LUS support and dims bins with n below 7, so a thin fit
is visible rather than silently shipped.

What to watch as samples are added:

- Does the LUS range with usable support widen past 10 to 14?
- Does `log_slope` stay near 0.72?
- Do the step factors stay near 0.24 and 0.73?
- Does any marker separate far enough from the pooled curve to justify a
  per-marker override?
