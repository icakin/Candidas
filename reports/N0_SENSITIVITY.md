# N0_SENSITIVITY — the N0 back-calculation, characterised, bounded and re-run

C2. Diagnostic plus a full downstream sensitivity analysis under three N0
treatments. **It changes no shipped default.** `N0_BACKPROJECT` is still `TRUE`.

Everything here was measured on this machine, from the C1 regeneration
(`runs/C1_reproduction/`) and the raw exports in `data/`. Tools and raw JSON:
`reports/tools/c2_*`. Figures: `reports/figures_N0/`.

---

## RECOMMENDATION, up front

**Keep ARM 1 (the shipped `N0 = N_inoc · e^(r·δ)`) for the paper as it stands,
and stop describing it as either "the artifact" or "clean". It is a
directionally correct correction with a measurable, quantified bias, and three
of the paper's four load-bearing claims survive it under every treatment
tested.** The fourth does not, and the paper should say so.

The reasoning, from the numbers rather than from convenience:

* **ARM 3 (`N0 = N_inoc`) is not a candidate and never was.** δ is a
  *detection* delay, not a lag — it is 58–82 min at every temperature with no
  trend (§1.2), whereas a biological lag would be ~3× shorter near 34 °C where
  growth is ~3× faster. The cells were present and growing through δ, so
  asserting zero growth is known to be wrong. It is reported throughout as a
  **lower bound on N0**, i.e. an **upper bound on per-cell respiration**.
* **ARM 2 (ramp-aware) is a bound, not an estimator**, because of the
  circularity (§2.4) — and it moves everything in the *same direction as ARM 1
  already goes*, not against it. It says the shipped correction is
  slightly too small below ~36 °C and slightly too large above it.
* So the two defensible treatments, ARM 1 and ARM 2, **bracket a narrow range**,
  and every claim that survives ARM 1 also survives ARM 2. ARM 3 sits well
  outside that bracket and is where things break.

**But one claim should be softened now, without waiting for new data.**
"*Clade IV is credibly the cheapest at fever*" rests entirely on the Clade III
vs Clade IV contrast, which resolves under ARM 1 (1.198 [1.031, 1.394]) and
ARM 2 (1.176 [1.029, 1.352]) but **not** under ARM 3 (1.100 [0.920, 1.319]).
That contrast is the thinnest in the set and it is the one the N0 treatment can
overturn. The species-level claim — all four *C. auris* clades pay credibly less
than *C. parapsilosis* — survives everywhere and is the safe headline.

**And one figure panel should be withdrawn or heavily caveated.** Fig 4e,
r(capacity, per-cell respiration) = −0.64, collapses to **−0.196 (p = 0.54)**
under ARM 3 (§3.5). Unlike everything else, this is not a matter of degree: the
correlation exists only while respiration carries the `e^(r·δ)` term, i.e. while
it is partly a function of growth rate. Capacity is fitted to growth. So Fig 4e
is at serious risk of being a correlation between two things that both contain
`r`. Fig 4d (r = −0.85, capacity vs fever tax) is stable across all three arms
and does not have this problem.

**What would actually settle it** is one endpoint cell count per vial (§5).
There is none in the repository (§4). Until then the estimand is not identified,
and the honest position is that ARM 1 and ARM 2 bracket it.

---

## 0. What δ is, and the identity that makes it matter

From the source, not assumed (`07_oxygen_fits.R:361, 922, 928, 935, 941`):

```
delta            = fit_start_time                      (the window start)
N0               = N_inoc * exp(r * delta)
C_tot_O2         = (K/r) * (exp(r*T_end) - 1)
biomass_integral = N0 * (exp(r*T_end) - 1) / r
R_per_cell       = C_tot_O2 / biomass_integral
```

The window length `T_end` **cancels exactly**, so

> **respiration = K / N0**, and therefore
> **log R = log K − log N_inoc − r·δ**

Verified numerically: the residual is constant at 30.745515174 with sd
**1.3 × 10⁻¹⁴** — exactly `log(MG_TO_FG · O2_TO_C_MASS · RQ · MIN_TO_H)`, the
unit conversion. Nothing else is in there. **Respiration is an explicit
decreasing function of the fitted growth rate.**

The same construction is in the published method
(`icakin/OxygenModel`, `CUE.R:222`), so everything below is **method-level, not
Candidas-specific**.

---

## 1. PART A — how big is the term, and what is δ?

### 1.1 Distributions (991 kept series, 6 groups, 12 temperatures)

| quantity | median | IQR | p95 | max |
|---|---:|---|---:|---:|
| δ (min) | 75.2 | 65.2 – 88.2 | 140.5 | 234.7 |
| r (min⁻¹) | 0.0088 | 0.0058 – 0.0120 | 0.0154 | 0.0208 |
| **r·δ** | **0.679** | 0.421 – 0.980 | 1.392 | 3.575 |
| **e^(r·δ)** | **1.97** | 1.52 – 2.67 | 4.02 | **35.7** |

So the median series has its N0 inflated ~2-fold, and hence its per-cell
respiration deflated ~2-fold, by this term alone.

**By temperature** — and this is the important shape:

| T (°C) | n | δ median | r·δ median | e^(r·δ) median | e^(r·δ) max |
|---:|---:|---:|---:|---:|---:|
| 22 | 89 | 58.1 | 0.226 | 1.25 | 2.27 |
| 26 | 90 | 74.2 | 0.541 | 1.72 | 7.02 |
| 30 | 90 | 80.2 | 0.832 | 2.30 | 7.79 |
| **32** | 90 | 78.2 | **0.978** | **2.66** | 14.46 |
| **34** | 89 | 76.2 | **0.980** | **2.66** | **35.68** |
| 36 | 87 | 76.2 | 0.938 | 2.55 | 16.79 |
| 40 | 78 | 79.2 | 0.743 | 2.10 | 7.69 |
| 44 | 47 | 82.2 | 0.422 | 1.52 | 2.80 |

The deflation **peaks at 32–34 °C and falls away at both ends**, tracking `r`.
It is not noise added to respiration; it is *the growth curve, inverted*,
imposed on respiration — deepest exactly where growth is fastest.

**At 22 / 34 / 44 °C, per taxon** (median e^(r·δ)):

| taxon | 22 °C | 34 °C | 44 °C |
|---|---:|---:|---:|
| Clade1 | 1.39 | 2.85 | 1.50 |
| Clade2 | 1.12 | 2.27 | 1.11 |
| Clade3 | 1.25 | 2.43 | 1.74 |
| Clade4 | 1.12 | **2.92** | 2.24 |
| glab | 1.49 | 3.58 | 1.33 |
| para | 1.32 | 2.66 | — (no 44 °C data) |

### 1.2 δ is instrumental, not biological

| test | result |
|---|---|
| δ ~ T | slope **+0.489 min/°C**, R² = **0.013**, p = 3.5e-4; Spearman ρ = +0.200 |
| δ ~ \|T − 20 °C\| | identical to the above — see the caveat |
| δ ~ r | R² = **0.002**, p = 0.15 (Spearman ρ = +0.144) |
| δ ~ taxon | ANOVA F(5,985) = 210, p = 1.5e-152 |
| δ ~ taxon, **excluding glabrata** | F(4,818) = **1.97, p = 0.097** |

Per-taxon mean δ: Clade1 69.9, Clade2 71.6, Clade3 69.8, Clade4 73.1,
para 73.2, **glab 125.6** (residual **+44.9 min**). The entire taxon effect is
glabrata, consistent with its documented signal-to-noise limit.

Median δ runs 58.1 min at 22 °C to 82.2 min at 44 °C — a 24 min spread across a
22 °C range, on a ~72 min baseline. **A biological lag would be ~3× shorter near
34 °C, where growth is ~3× faster. It is not; if anything it is slightly
longer.**

> **CAVEAT on the \|T − T_ambient\| regression.** Every temperature in the design
> is above the assumed 20 °C ambient, so `|T − 20|` is a strictly increasing
> function of `T` and the two regressions are **numerically identical**. This
> design cannot separate them. The instrumental reading is therefore tested
> properly in PART B, on the transient's own amplitude and time constant, where
> it is not confounded.

**δ is not the automatic O2 peak.** It equals `peak_time` in only **16 of 991**
series (1.6%); the median difference is +2.0 min but the range is
−108 to +148 min. δ is the *hand-chosen fit-window start* from
`04_trim_selector.R`. Both are "time from liquid-in to the start of the fit",
so the argument is unaffected, but the provenance should be stated correctly.

### 1.3 How much of E_R is the correction?

Removing the term makes E_R **identical to E_K**, the N0-free volumetric slope,
by construction (`log R_nobp = log K − log N_inoc`, and `N_inoc` has no
temperature dependence). Full range, 22–44 °C:

| taxon | E_K = E_R(removed) | E_R shipped | artefact (E_K − E_R) | median r·δ |
|---|---:|---:|---:|---:|
| Clade1 | 0.625 ± 0.030 | 0.518 ± 0.033 | **0.107** | 0.675 |
| Clade2 | 0.545 ± 0.030 | 0.416 ± 0.035 | **0.129** | 0.429 |
| Clade3 | 0.592 ± 0.020 | 0.404 ± 0.029 | **0.189** | 0.733 |
| Clade4 | 0.597 ± 0.021 | 0.291 ± 0.032 | **0.306** | 0.893 |
| glab | 0.825 ± 0.040 | 0.694 ± 0.080 | 0.131 | 0.797 |
| para | 0.547 ± 0.043 | 0.339 ± 0.060 | **0.207** | 0.656 |

The artefact is **ordered by growth rate** — largest for Clade4, which has the
largest median r·δ. This is the inverted-growth-curve mechanism, visible
directly.

**The between-taxon ordering does not survive removal.**

* shipped: Clade4 < para < Clade3 < Clade2 < Clade1 < glab
* removed: Clade2 < para < Clade3 < Clade4 < Clade1 < glab
* Spearman between the two orderings: **+0.486 (p = 0.33)**
* between-taxon spread of E_R: **0.403 eV → 0.280 eV**
* **Clade IV moves from LOWEST E_R (0.291) to 4th of 6 (0.597).**

On the **rising limb (T ≤ 34 °C)** it is starker still: shipped E_R spans
0.054–0.358 eV, and Clade4's is **0.054 ± 0.057 — indistinguishable from zero**.
With the term removed all six lie in 0.55–0.86.

*Framing, as the brief requires:* this measures **how much of E_R is the
correction**. It is not evidence the correction is wrong. §1.2 shows the cells
really were growing through δ.

**Validation.** These cheap OLS Boltzmann–Arrhenius slopes track C1's
hierarchical brms E_R to within **0.024 eV** for every taxon (Clade1 0.5176 vs
0.5193; Clade4 0.2906 vs 0.3006), so PART A's arithmetic is a faithful stand-in
for the fitted model.

---

## 2. PART B — the transient, and the two errors it creates

### 2.1 The transient itself

`O2(t) = O2_eq − A · exp(−t/τ)` fitted to every pre-peak segment:
**1074 of 1080 converge** (6 do not rise), median R² **0.9988**
(IQR 0.9980–0.9994), median RMSE 0.027 mg/L.

| T (°C) | \|T−20\| | A (mg/L) median [IQR] | τ (min) median [IQR] | O2_eq (mg/L) | settled at the peak |
|---:|---:|---|---|---:|---:|
| 22 | 2 | 1.04 [0.54, 1.43] | 17.0 [14.1, 19.3] | 10.74 | 98.0 % |
| 28 | 8 | 2.78 [2.25, 3.27] | 18.4 [15.2, 22.8] | 10.45 | 98.3 % |
| 34 | 14 | 4.01 [3.39, 4.47] | 20.0 [16.7, 23.0] | 9.02 | 97.7 % |
| 38 | 18 | 4.98 [4.26, 5.27] | 20.1 [18.3, 22.5] | 8.96 | 97.6 % |
| 44 | 24 | 5.70 [5.21, 6.08] | 20.0 [17.4, 22.5] | 8.97 | 98.1 % |

**It is physico-chemical.** The amplitude scales with the temperature step
(`A ~ |T−20|`: slope +0.2035 mg/L/°C, **R² = 0.822**, Spearman ρ = +0.910), and
against the Benson–Krause solubility curve,
**Pearson r(A, |sol(T) − sol(20)|) = +0.996, p = 9.5 × 10⁻¹²**, with the ratio
settling at ~1.9–2.0 for T ≥ 32 °C.

**τ does not differ by taxon.** ANOVA on per-temperature residuals excluding
glabrata: **F(4,889) = 1.99, p = 0.094**. glabrata is +2.8 min on a ~19 min
baseline. Nothing here is biological.

**O2_eq vs solubility.** Fitted O2_eq falls with T (slope −0.091 mg/L/°C)
essentially in parallel with solubility (−0.121); **Pearson r = +0.892**
across the 12 temperatures. But it sits **1.20–1.48× above air-saturation**,
with the excess growing with T.

**Two refinements to the brief's premise, both worth recording.**

1. **τ does *not* scale with distance from ambient** — it is 17.0–20.4 min at
   every temperature (slope +0.108 min/°C, R² = 0.027). What lengthens at high
   temperature is the time to the *peak*, and only because a larger rise takes
   longer to be overtaken by the biological drawdown. The 37 / 51 / 70 min
   figures in the brief are peak times, not time constants.
2. **The sign is inverted.** Warming should *lower* dissolved O2; the record
   rises. Combined with O2_eq sitting above air-saturation, this is the
   **optode's own temperature compensation settling**, not gas
   re-equilibration — and §2.3 identifies the mechanism exactly. Either way the
   driver is the temperature step, not the cells.

### 2.2 (B1) Fit-window contamination — real, temperature-dependent, and the wrong sign to explain anything

**How it was measured.** The fitted transient is extrapolated forward,
subtracted from the observed trace over the actual fit window, and 07's own
oxygen model `O2 = O2_0 + (K/r)(1 − e^(rt))` is **refitted** to the
decontaminated trace with 07's own bounds. The refit first reproduces the
**shipped** K from the raw window (median |rel diff| 4.9e-3; 66 % within 1 %),
so the comparison is like for like.

> An earlier pass compared the residual abiotic *slope* at the window start to K
> and got 60–79 %. That is wrong — it sets an instantaneous rate against a
> whole-window parameter — and is not what is reported.

**Sign:** the transient *rises*, so residual drift inside the window masks part
of the drawdown and makes K **too small**.

| T (°C) | settled at window start | abiotic drift (mg/L) | as % of drawdown | **K under-read**, τ q25 → **fitted τ** → τ q75 |
|---:|---:|---:|---:|---|
| 22 | 97.82 % | 0.016 | 0.52 % | 0.21 % → **0.54 %** → 0.99 % |
| 28 | 98.32 % | 0.047 | 1.02 % | 0.76 % → **2.16 %** → 5.34 % |
| 34 | 97.97 % | 0.075 | 1.83 % | 2.42 % → **5.02 %** → 10.46 % |
| 36 | 97.82 % | 0.093 | 2.40 % | 3.02 % → **6.41 %** → 13.48 % |
| 40 | 97.91 % | 0.106 | 2.46 % | 3.65 % → **6.63 %** → 11.09 % |
| 44 | 98.17 % | 0.097 | 2.27 % | 2.14 % → **4.01 %** → 7.35 % |

It **does** grow with temperature (+0.288 %-points/°C, p = 1.3e-23,
Spearman ρ = +0.521), peaking near 40 °C.

**But it cannot account for the monotonic rise of respiration to 44 °C, because
correcting it makes the rise steeper, not shallower.** E_K after correction, at
the fitted τ (and the whole τ-IQR band keeps this sign):

| taxon | E_K raw | at τ q25 | at fitted τ | at τ q75 |
|---|---:|---:|---:|---:|
| Clade1 | 0.621 | 0.657 | **0.656** | 0.691 |
| Clade2 | 0.546 | 0.565 | **0.571** | 0.591 |
| Clade3 | 0.584 | 0.597 | **0.601** | 0.614 |
| Clade4 | 0.569 | 0.583 | **0.591** | 0.603 |
| glab | 0.847 | 0.853 | **0.851** | 0.858 |
| para | 0.540 | 0.558 | **0.574** | 0.591 |

This is a **model-based sensitivity band, not a point estimate** — there are no
cell-free wells anywhere in the design (row D and column 6 are "No Sensor"), so
the transient inside the window is extrapolated, never observed.

### 2.3 The mechanism, from the raw exports

While searching `data/` for endpoint biomass (§4) the raw PreSens exports turned
out to carry two temperature columns that `01_convert_xlsx.R` discards:

* **`Tm [°C]`** — **constant for the whole run**, and equal to the file-name set
  point in **72 of 72 plates**. It is the value the operator typed into the SDR
  software, and (with `Temperature comp. : ON` in the calibration header) it is
  what the optode's temperature compensation uses.
* **`T_internal [°C]`** — a **real measured temperature trace**.

So the compensation assumed the set point for the whole run while the plate was
measurably 1.6–2.9 K cooler for the first hour. That is a complete account of
the pre-peak transient, of why its amplitude scales with |T − ambient|, and of
why O2_eq sits above air-saturation. The calibration header also records
`O2_T0 = O2_T2nd = 23.2 °C` — the sensor was two-point calibrated at 23.2 °C.

### 2.4 (B2) Thermal-ramp error — measured, not assumed

`reports/tools/c2_extract_Tinternal.R` recovers `T_internal` from all 72 plates
(99,120 rows), which turns the thermal ramp from an assumption into a
measurement.

| set point | t=0 | t=20 | t=40 | t=75 | t=120 | settled | mean deficit over 0–75 min |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 22 | 22.24 | 22.13 | 22.43 | 22.59 | 22.63 | 22.67 | **0.30 K** |
| 34 | 32.53 | 31.43 | 32.64 | 33.48 | 33.64 | 33.75 | **1.60 K** |
| 44 | 41.97 | 39.32 | 41.45 | 43.01 | 43.39 | 43.57 | **2.86 K** |

* the plate **does not start at ambient** — at 44 °C it begins 1.6 K below its
  own asymptote, not 24 K below;
* it **dips** as the load equilibrates and recovers over roughly an hour;
* at 22 °C there is effectively no ramp at all.

Ramp-aware N0 = `N_inoc · exp(∫₀^δ r(T(t)) dt)`, with `T(t)` the measured trace
and `r(T)` each isolate's own Sharpe–Schoolfield TPC scaled through that series'
own fitted `r` (so only the *shape* of the TPC is used).

| T (°C) | mean T over 0..δ | N0_ramp / N0_current, median [IQR] | r·δ current → ramp |
|---:|---:|---|---|
| 22 | 22.37 | 0.9920 [0.9898, 0.9943] | 0.226 → 0.220 |
| 30 | 28.37 | 0.9208 [0.9143, 0.9321] | 0.832 → 0.761 |
| **32** | 30.00 | **0.9150** [0.8951, 0.9315] | 0.978 → 0.884 |
| 36 | 34.27 | 0.9802 [0.9585, 1.0232] | 0.938 → 0.973 |
| 40 | 38.22 | 1.1394 [1.0918, 1.2022] | 0.743 → 0.866 |
| 44 | 41.44 | **1.1996** [1.1395, 1.2724] | 0.422 → 0.665 |

**The ratio is not monotone, and the optimum is why.** Below ~36 °C the ramp
sits below the set point, so the cells grew *slower* than `r` and the current N0
**over-corrects**. Above it, the ramp passes *through* the optimum, so the cells
grew *faster* than at the supra-optimal set point and the current N0
**under-corrects**. Net effect on E_R: **−0.009 to −0.102 eV** in every taxon.

*The modelled version this replaces* (start at 20 °C ambient, τ = 19 min)
**overstated the ramp roughly twofold**: 0.841 vs 0.941 at 34 °C, 0.885 vs 0.980
at 36 °C. Kept for comparison as `c2_ramp_n0_factors_modelled.csv`.

**The circularity, handled explicitly.** The growth TPC is **N0-invariant** —
`r` comes from the oxygen curve shape and N0 enters only respiration = K/N0 — so
there is no feedback through the growth side. The real coupling is the **342 of
782 series (43.7 %)** still more than 0.25 K below settled when their fit window
opens, where the fitted `r` is itself a ramp average. One iteration was done;
the second-iteration change in the ratio is median 0.0048, p95 0.065, max 0.360.
**ARM 2 is therefore a BOUND on the size and direction of the ramp error, not a
new estimator.**

> **Caveat.** `T_internal` is the reader's internal sensor, not a probe in the
> liquid, so it remains a *proxy* for vial temperature. It is a **measured**
> proxy, which the previous assumption was not.

---

## 3. PARTs C & D — three arms, end to end

| | treatment | tree | wall time |
|---|---|---|---:|
| **ARM 1** | `N0 = N_inoc · e^(r·δ)` (shipped) | `runs/C2_arm_current/` | 10 m 33 s |
| **ARM 2** | `N0 = N_inoc · e^(r·δ) · f_ramp` | `runs/C2_arm_ramp/` | 10 m 23 s |
| **ARM 3** | `N0 = N_inoc` — **LOWER BOUND on N0, not a candidate** | `runs/C2_arm_nobp/` | 11 m 07 s |

Each is a complete 07 → 08 → 09 → 10 → 12 → 13 with a **full brms refit**.
Held fixed across arms: the committed trim windows, exclusions, inoculum and
cell-carbon constants; every brms formula, prior and the seed 1234; and the
etc-GEM capacity values (`runs/C1/supp_data/` for every arm — C1 established
the etc-GEM is not reproducible, so letting capacity vary would confound the
comparison). 02/03/06/11 are N0-independent and their outputs are carried over
from `runs/C1_reproduction/`, so all three arms literally share the same inputs.

brms per arm: `iter 4000, warmup 1000, chains 4, seed 1234, adapt_delta 0.99 /
0.95, max_treedepth 12`. **max R̂ = 1.00297, min ESS bulk = 2601, no divergent
transitions**, identically in all three arms (the growth fit is bit-identical
across arms, so these are the same numbers, not merely similar ones).

### 3.0 The reproducibility floor (PART C0) — established first

`runs/C2_arm_current` vs `runs/C1_reproduction`: same treatment, different run.

**35 of 36 comparable tables are IDENTICAL (max |diff| < 1e-10)** — including
`fig_values.csv`, `fig_contrasts.csv`, `carbon_tax_group.csv` and
`carbon_tax_isolate.csv`. The single exception is `carbon_tax_curves.csv` at
**0.165 absolute / 2.11 % relative**, inside C1's measured *unseeded* floor of
0.205 / 2.21 %.

**Seed: `CANDIDAS_SEED = 20260726`**, applied identically in all three arms, at
`10_carbon_tax.R:217` and `12_main_figures.R:358, 621, 674, 1152`. Seeding
**moved no point estimate**: every quantity in the tables below has a floor of
**exactly zero**, and only the Fig 2a curves have a non-zero floor (2.1 %).
Anything smaller than that is not evidence; nothing reported below is that
small.

### 3.1 INVARIANCE CHECK — PASS, exactly

N0 touches only respiration, so growth must not move. Maximum |difference|
between every pair of arms:

| quantity | ramp vs current | nobp vs current |
|---|---:|---:|
| `r` | **0.000e+00** | **0.000e+00** |
| `K` | **0.000e+00** | **0.000e+00** |
| `growth_fgC_h` | **0.000e+00** | **0.000e+00** |
| `growth_C_per_C_h` | **0.000e+00** | **0.000e+00** |
| growth Sharpe-Schoolfield posterior summary | **0.000e+00** | **0.000e+00** |
| per-isolate SS growth coefficients | **0.000e+00** | **0.000e+00** |
| Arrhenius growth coefficients | **0.000e+00** | **0.000e+00** |
| E_G | **0.000e+00** | **0.000e+00** |
| growth T_opt | **0.000e+00** | **0.000e+00** |

**The switch does not leak.** (C1 saw its growth fit move between its own two
runs because the *input* `se_y` differed in the last bits across a version
change; here the inputs are literally the same files, and the result is exact
zero — so that effect is absent, not merely small.)

### 3.2 Per-taxon quantities

**E_R (eV)** — the quantity the whole question is about:

| taxon | ARM 2 ramp | ARM 1 current | ARM 3 no-backproj |
|---|---|---|---|
| Clade1 | 0.442 [0.372, 0.514] | 0.519 [0.443, 0.594] | 0.654 [0.591, 0.716] |
| Clade2 | 0.352 [0.273, 0.433] | 0.414 [0.331, 0.497] | 0.553 [0.484, 0.620] |
| Clade3 | 0.329 [0.259, 0.400] | 0.401 [0.326, 0.479] | 0.600 [0.541, 0.660] |
| Clade4 | 0.238 [0.167, 0.312] | 0.301 [0.226, 0.380] | 0.622 [0.561, 0.682] |
| para | 0.283 [0.188, 0.377] | 0.363 [0.264, 0.463] | 0.551 [0.478, 0.628] |

**E_G (eV)** — identical in all three arms, as it must be: Clade1 0.914
[0.795, 1.055], Clade2 1.143 [1.010, 1.295], Clade3 0.947 [0.832, 1.074],
Clade4 1.019 [0.900, 1.154], para 0.830 [0.699, 0.985].

**ΔE = E_G − E_R (eV)**, and whether it excludes zero:

| taxon | ARM 2 ramp | ARM 1 current | ARM 3 no-backproj |
|---|---|---|---|
| Clade1 | 0.472 [0.331, 0.631] ✔ | 0.395 [0.254, 0.556] ✔ | 0.261 [0.125, 0.414] ✔ |
| Clade2 | 0.791 [0.634, 0.962] ✔ | 0.730 [0.570, 0.900] ✔ | 0.590 [0.441, 0.758] ✔ |
| Clade3 | 0.619 [0.482, 0.763] ✔ | 0.547 [0.407, 0.695] ✔ | 0.347 [0.219, 0.485] ✔ |
| Clade4 | 0.781 [0.641, 0.934] ✔ | 0.718 [0.576, 0.872] ✔ | 0.397 [0.265, 0.548] ✔ |
| para | 0.547 [0.382, 0.727] ✔ | 0.468 [0.303, 0.653] ✔ | 0.279 [0.127, 0.451] ✔ |

**ΔE excludes zero in every taxon in every arm.**

**T_opt(CUE) (°C)**, and P(T_opt < 37 °C):

| taxon | ARM 2 ramp | ARM 1 current | ARM 3 no-backproj |
|---|---|---|---|
| Clade1 | 32.01 [30.90, 33.14] | 31.27 [30.03, 32.45] | 29.63 [28.01, 30.90] |
| Clade2 | 33.25 [32.38, 34.15] | 32.98 [32.08, 33.88] | 32.30 [31.42, 33.20] |
| Clade3 | 33.35 [32.32, 34.35] | 32.77 [31.72, 33.83] | 30.88 [29.72, 31.96] |
| Clade4 | 34.58 [33.61, 35.57] | 34.11 [33.08, 35.13] | 31.21 [30.10, 32.35] |
| para | 31.78 [30.73, 32.81] | 31.30 [30.20, 32.40] | 29.84 [28.44, 31.01] |

**P(T_opt(CUE) < 37 °C) = 1.0000 for every taxon in every arm.**

**Carbon tax:**

| taxon | tax37 ramp / current / nobp | tax40 ramp / current / nobp |
|---|---|---|
| Clade1 | 1.233 / 1.298 / 1.446 | 1.756 / 1.899 / 2.221 |
| Clade2 | 1.264 / 1.301 / 1.401 | 2.136 / 2.247 / 2.542 |
| Clade3 | 1.130 / 1.171 / 1.327 | 1.512 / 1.608 / 1.959 |
| Clade4 | 1.052 / 1.073 / 1.269 | 1.285 / 1.341 / 1.777 |
| para | 1.512 / 1.597 / 1.850 | 2.875 / 3.127 / 3.879 |

**CUE by clade** (median, peak):

| taxon | ARM 1 current | ARM 2 ramp | ARM 3 no-backproj |
|---|---|---|---|
| Clade1 | 0.340, 0.379 | 0.335, 0.379 | **0.220, 0.241** |
| Clade2 | 0.211, 0.275 | 0.210, 0.278 | **0.151, 0.194** |
| Clade3 | 0.411, 0.463 | 0.411, 0.468 | **0.269, 0.300** |
| Clade4 | 0.486, 0.543 | 0.488, 0.547 | **0.296, 0.330** |
| para | 0.321, 0.380 | 0.319, 0.384 | **0.209, 0.238** |

Absolute CUE is not identifiable anyway (it scales with cell volume and
inoculum), which is why the paper's claims live in ratios — but ARM 3 lowers it
by ~35 % across the board.

### 3.3 The 10 pairwise fever contrasts

`*` = 95 % CrI excludes 1.

| a | b | ARM 1 current | ARM 2 ramp | ARM 3 no-backproj |
|---|---|---|---|---|
| Clade4 | para | 0.430 [0.326, 0.551] * | 0.447 [0.341, 0.565] * | 0.459 [0.346, 0.595] * |
| Clade3 | para | 0.514 [0.386, 0.668] * | 0.526 [0.398, 0.673] * | 0.504 [0.376, 0.657] * |
| Clade1 | para | 0.607 [0.452, 0.797] * | 0.611 [0.460, 0.789] * | 0.573 [0.426, 0.756] * |
| Clade2 | para | 0.718 [0.523, 0.956] * | 0.741 [0.547, 0.973] * | 0.656 [0.479, 0.871] * |
| Clade1 | Clade2 | 0.846 [0.683, 1.043] | 0.823 [0.672, 1.009] | 0.873 [0.698, 1.091] |
| Clade1 | Clade3 | 1.181 [0.983, 1.414] | 1.161 [0.980, 1.369] | 1.134 [0.927, 1.386] |
| **Clade3** | **Clade4** | **1.198 [1.031, 1.394] \*** | **1.176 [1.029, 1.352] \*** | **1.100 [0.920, 1.319]** |
| Clade2 | Clade3 | 1.395 [1.143, 1.711] * | 1.410 [1.172, 1.704] * | 1.298 [1.050, 1.615] * |
| Clade1 | Clade4 | 1.416 [1.200, 1.666] * | 1.365 [1.172, 1.592] * | 1.247 [1.028, 1.519] * |
| Clade2 | Clade4 | 1.674 [1.390, 2.028] * | 1.662 [1.383, 1.987] * | 1.429 [1.162, 1.760] * |

**Resolved: ARM 1 8/10 · ARM 2 8/10 · ARM 3 7/10.**
The one that flips is **Clade III vs Clade IV**.

### 3.4 Spearman ρ (environmental optimum vs fever cost)

| | ρ | 95 % CrI | P(ρ < 0) |
|---|---:|---|---:|
| ARM 1 current | −1.00 | [−1.00, −0.70] | 1.000 |
| ARM 2 ramp | −1.00 | [−1.00, −0.70] | 1.000 |
| ARM 3 no-backproj | −1.00 | [−1.00, −0.70] | 1.000 |
| *baseline: the committed brms fits (C1, twice)* | −1.00 | [−1.00, −0.70] | 1.000 |
| **the manuscript states** | **−0.90** | **[−1.00, −0.30]** | — |

**Invariant to the N0 treatment.** The manuscript/result mismatch is a C1
finding, restated here so it is not lost; **no fourth number is introduced**.

### 3.5 Fig 4d and 4e — capacity correlations

Capacity is held fixed at `runs/C1/supp_data/` in all three arms, so only the
economics axis moves.

| arm | **4d** r(capacity, fever tax at 40 °C) | **4e** r(capacity, per-cell respiration) |
|---|---|---|
| ARM 1 current | −0.847 (p = 5.0e-4) | −0.620 (p = 3.2e-2) |
| ARM 2 ramp | −0.874 (p = 2.1e-4) | −0.603 (p = 3.8e-2) |
| **ARM 3 no-backproj** | **−0.819 (p = 1.1e-3)** | **−0.196 (p = 0.54)** |

**4d is stable. 4e collapses.** Per-cell respiration under ARM 3 is `K/N_inoc`,
with no `r` in it; capacity is fitted to growth. The correlation exists only
while respiration carries the `e^(r·δ)` term. This is the single most
treatment-sensitive number in the paper.

> Reported for **arm comparison only**. The capacity axis itself is unresolved
> pending the etc-GEM re-run — C1 established that the committed etc-GEM outputs
> do not reproduce (anchor change; see `reports/REPRODUCTION.md` §4.2).

### 3.6 Verdict on each named claim

| claim | ARM 1 current | ARM 2 ramp | ARM 3 no-backproj (lower bound) |
|---|---|---|---|
| **"growth is ~2× more temperature-sensitive than respiration"** | **SURVIVES** (E_G/E_R 1.76–3.39, median 2.36) | **SURVIVES** (2.07–4.28, median 2.94) | **SURVIVES as a direction, WEAKENS as "2×"** (1.40–2.07, median 1.58; ΔE still excludes zero in all five taxa) |
| **"Clade IV is credibly the cheapest at fever"** | **SURVIVES** (III vs IV 1.198 [1.031, 1.394]) | **SURVIVES** (1.176 [1.029, 1.352]) | **FAILS** (1.100 [0.920, 1.319] — CrI spans 1) |
| **"all four *C. auris* clades pay credibly less than *C. parapsilosis*"** | **SURVIVES** (4/4) | **SURVIVES** (4/4) | **SURVIVES** (4/4) |
| **"every economic optimum lies below 37 °C"** | **SURVIVES** (P = 1.0000 ×5) | **SURVIVES** (P = 1.0000 ×5) | **SURVIVES** (P = 1.0000 ×5) |

Every difference above is far larger than the measured floor of exactly zero
(§3.0). The only quantity with a non-zero floor is the Fig 2a curve band
(2.1 %), and no claim rests on it.

---

## 4. PART E — what would settle it, and what the repository has

### 4.1 The two available back-calculations

|  | **FORWARD (current)** | **BACKWARD (proposed)** |
|---|---|---|
| formula | `N0 = N_inoc · exp(r·δ)` | `N0 = N_end · exp(−r·T_window)` |
| assumes | growth at rate `r` through the equilibration and thermal ramp — an interval the oxygen record **cannot see**, and during which the cells were measurably **not at temperature** (§2.4: 1.6–2.9 K below settled for the first hour) | exponential growth **over the fit window only** |
| extra assumption over the model | **yes** — extrapolates the growth model backwards, outside the fitted window | **none.** The oxygen model `O2 = O2_0 + (K/r)(1 − e^(rt))` *already* assumes exponential growth across exactly that window |
| what it needs | nothing extra | **one endpoint count per vial**, taken **after** the run |
| disturbs the sealed measurement? | no | **no** — the count is taken after the seal is broken at the end |

The backward route **adds no assumption the model does not already make**. That
is its whole appeal: it replaces an extrapolation into an unobserved interval
with a measurement at a moment the experiment is already over.

**Running BOTH yields the missing quantity empirically.** The ratio

```
    [N_end · exp(−r·T_window)]  /  [N_inoc · exp(r·δ)]
```

*measures* how much growth actually occurred during equilibration — the number
this entire report has had to bound rather than observe. It would settle the
question in one experiment, and it would settle it for the published method too,
not just for this dataset.

> **Not proposed, and deliberately so:** measuring cell density in the vial at
> the window start. It cannot be done without breaking the seal and destroying
> the closed-system oxygen measurement. The back-calculation is a considered
> design choice; the question is only *which* one is best supported.

### 4.2 What endpoint biomass exists in the repository

**None.** Searched exhaustively:

| searched | result |
|---|---|
| `data/*_Oxygen.csv` — every distinct column name | `Time`, `T`, `OTU<n>_R<n>` only. No biomass column. |
| `data/*.xlsx` — all 72 raw exports, all 32 columns | `Date/Time`, `Time/Min.`, wells `A1`–`D6`, `Tm [°C]`, `p [mbar]`, `Salinity [g/1000g]`, `T_internal [°C]`, `Error`. **No optical density, no cell count, no biomass.** |
| `data/` subdirectories | only `data/expression/` (GEO GSE165762 RNA-seq — not biomass) |
| repo-wide grep for haemocytometer / flow cytometry / cell count / final OD / OD600 | hits only in `06_inoculation.R` and `07_oxygen_fits.R` *comments*, and in the C2 prompts themselves. No data. |
| `results/tables/otu_inoc.csv` | **INOCULUM only** — `N_inoc_cells_per_L` per isolate, i.e. the *start*, which is what the forward route already uses. Nothing at the end. |

So the backward route **cannot be run retrospectively on this dataset**. The
counts were never taken.

> **Flagged, not attempted:** the OxygenModel paper (`icakin/OxygenModel`), which
> shares this exact construction at `CUE.R:222`, is reported to have
> flow-cytometry data for its 15 bacterial taxa. If those are *endpoint* counts
> per vial, that dataset could support the both-routes test retrospectively and
> would resolve the method-level question without new bench work. **Worth
> asking for. Not obtained here, and no claim is made about what it would show.**

### 4.3 The two assay additions, concretely

**(1) Cell-free medium-only wells at EVERY temperature.** Measure the transient
instead of extrapolating it — this is the one thing that would turn §2.2 from a
sensitivity band into a measurement.

* *What:* wells containing medium only, no inoculum, otherwise identical
  (same volume, same seal, same optode, same run).
* *How many:* **2 per plate**, i.e. 2 × 12 temperatures × 6 groups = 24 wells,
  or 2 × 12 = 24 if one group's plates suffice. The plate already has spare
  capacity — **row D and column 6 are currently "No Sensor"**, 9 unused
  positions per 24-well plate. This costs no extra plates and no extra runs.
* *When:* every run, from now on.
* *What it buys:* the abiotic O2 trajectory observed directly through the entire
  fit window, so K can be corrected rather than bounded, and the τ-IQR band in
  §2.2 collapses to a measurement.

**(2) One endpoint count per vial.** Anchor N0 backwards.

* *What:* a cell count (haemocytometer or flow cytometry) from each vial
  immediately after the seal is broken at the end of the run.
* *How many:* **one per vial — all 15 per plate**, since the whole point is the
  per-series pairing with that vial's own `r` and window.
* *When:* at the end of every run, after the oxygen record is complete. **The
  sealed measurement is never disturbed.**
* *What it buys:* `N0 = N_end · exp(−r·T_window)`, and — run alongside the
  forward route — the empirical ratio in §4.1 that measures growth during
  equilibration.

Also, at no cost: **keep `T_internal` and `Tm`.** `01_convert_xlsx.R` discards
both. `T_internal` is what made §2.4 a measurement rather than an assumption,
and it should be carried into `Oxygen_All_Long.csv` alongside `Oxygen`.

### 4.4 Recommendation

**Build the paper on ARM 1, the shipped treatment, and say what it costs.**

* ARM 1 and ARM 2 bracket the defensible range and agree on every claim.
  ARM 2 says the shipped correction is a little too strong below ~36 °C and a
  little too weak above it, with E_R shifts of −0.009 to −0.102 eV — real, but
  not claim-changing.
* ARM 3 is a lower bound, not an alternative, and is reported as such. It is
  nonetheless the useful stress test, and two things break under it:
  **the Clade III vs Clade IV contrast**, and **Fig 4e**.
* Therefore: soften "Clade IV is credibly the cheapest at fever" to the
  species-level claim that survives everywhere, and either withdraw Fig 4e or
  state plainly in its legend that the correlation depends on the N0 treatment
  (−0.62 under the shipped one, −0.20 and non-significant under the lower
  bound).

**Is the estimand adequately identified? No.** Not by any of the three arms.
δ is real, growth through it is real, and its magnitude is *bounded* by this
analysis but not *measured*. The fix is cheap, uses plate capacity that already
exists, and is specified in §4.3. Until then, quote ARM 1 with the ARM 2 bracket
and treat the Clade III/IV ordering and Fig 4e as provisional.

**A note on scope.** `N0_BACKPROJECT` has **not** been changed, and neither has
anything else on the constraint list. The only source edits C2 made are the arm
plumbing (`CANDIDAS_N0_MODE`, which reproduces the shipped pipeline exactly when
unset), the seed helper, and the two comment blocks in §5.

---

## 5. The contradiction, adjudicated

Two comment blocks in this repository flatly contradicted each other on this
switch, and **neither cited an analysis**.

**`scripts/config.R`, at `N0_BACKPROJECT`, verbatim:**

> "TRUE = N0 = N_inoc * exp(r * delta) … Back-projects biomass across the
> (largely lag-phase) interval before the fit window; **the exp(r * delta) term
> injects the inverted growth curve and fit noise into respiration -> no thermal
> optimum, noisy.**
> FALSE = N0 = N_inoc (delta = 0) … **This removes the artifact so respiration
> shows a clean TPC comparable to growth.**"

**`scripts/04_trim_selector.R`, at `STAB_MIN_START_SLOPE_FRAC`, verbatim:**

> "CORRECTION (this comment used to blame the N0 back-projection and recommend
> N0_BACKPROJECT <- FALSE; that was wrong). glab's bad K is a windowing problem
> after all - but it is the window END, not the start. … **Do NOT set
> N0_BACKPROJECT <- FALSE.**"

**Which is supported?** Partly each, and neither as written.

| assertion | verdict | evidence |
|---|---|---|
| config.R: "injects the inverted growth curve … into respiration" | **SUPPORTED, and stronger than stated** | `log R = log K − log N_inoc − r·δ` exactly (§0). e^(r·δ) peaks at 32–34 °C — the growth optimum — at 2.66× and reaches 35.7× (§1.1). The E_R artefact is ordered by growth rate across taxa (§1.3). |
| config.R: "and fit noise" | **SUPPORTED** | δ is the hand-chosen window start, and `r` is a fitted parameter; both enter multiplicatively. |
| config.R: "no thermal optimum, noisy" | **NOT TESTED HERE** | This is a claim about the *shape* of the respiration TPC, not about E_R. LOO already prefers Arrhenius over Sharpe-Schoolfield for respiration; C2 did not re-open that. |
| config.R: FALSE "removes the artifact" | **NOT SUPPORTED** | δ is a detection delay, not a lag (§1.2). FALSE asserts zero growth over an interval where growth demonstrably occurred; it is a **lower bound on N0**, not a corrected estimator. |
| config.R: FALSE gives "a clean TPC" | **MISLEADING** | It gives E_R ≡ E_K, the raw volumetric slope. "Clean" only in the sense of *containing no correction at all*. |
| 04: glabrata's bad K is a window-**END** problem, fixed by capping rt | **SUPPORTED, and untouched by C2** | Not what C2 tested. The rt diagnosis stands on its own evidence (sd(log K) 0.58 → 0.37). |
| 04: "Do NOT set N0_BACKPROJECT <- FALSE" | **WITHDRAWN as written** | It is an instruction, not a finding. Fixing glabrata's windows says nothing about the other five taxa. The narrower true statement — FALSE is a lower bound, so it should not be *adopted* — is a reason not to adopt it, **not** a reason not to *test* it. |

**Both blocks have been rewritten to point at this report.** The edits are
**comments only**; `git diff` on the two files contains no non-comment line, and
`N0_BACKPROJECT` is still `TRUE`.

---

## 6. Reproducing this

```bash
# PART A/B diagnostics (read-only)
cauris_etcgem/.venv/bin/python reports/tools/c2_partA_quantify_delta.py
cauris_etcgem/.venv/bin/python reports/tools/c2_partB_transient.py
Rscript                        reports/tools/c2_extract_Tinternal.R
cauris_etcgem/.venv/bin/python reports/tools/c2_partB5_ramp_measured.py

# PART C: the three arms  (~32 min total)
bash scripts/run_c2_arms.sh                       # or: bash scripts/run_c2_arms.sh ramp

# PART D: comparison + figures
cauris_etcgem/.venv/bin/python reports/tools/c2_compare_arms.py
cauris_etcgem/.venv/bin/python reports/tools/c2_figures.py
```

| output | what |
|---|---|
| `reports/figures_N0/c2_arm_comparison.png` | E_R, ΔE, T_opt(CUE), tax40 — all taxa, all arms, with CrIs |
| `reports/figures_N0/c2_transient_diagnostics.png` | A vs T, τ vs T, and the measured thermal ramp |
| `reports/figures_N0/FIG{1,2}_*_<arm>.png`, `FIG_MODEL_<arm>.png` | each arm's own Fig 1, Fig 2 and Fig 4d/4e, taken from that arm's output tree (never from `results/figures/`, which C1 found holds 7 stale duplicates) |
| `reports/tools/c2_*.json`, `c2_partB_per_series.csv`, `c2_Tinternal.csv` | every number above |
