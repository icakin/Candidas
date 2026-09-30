# R1 — final verdict: not salvageable from these data (and a warning for Figs 1–2)

## What R1 asked
Constrain the etcGEM with measured oxygen and predict growth. This requires oxygen and
growth to be **independent measurements**. They are not.

## The exact structure of the data (verified, floating-point exact)
Each oxygen trace is fitted with ONE 3-parameter model (`07_oxygen_fits.R`, nlsLM on raw
oxygen):

    O2(t) = O2_0 + (K/r) * (1 - exp(r*t))          parameters: r, K, O2_0

The script states the premise plainly: *"The BEND in the O2 curve IS the growth."*
From that single fit:

- `growth_C_per_C_h  ==  r * 60`            (identity; max rel. deviation 7×10⁻¹⁴)
- `biomass_integral  ==  N0*(exp(r*T_end)-1)/r`   (identity; 5×10⁻¹⁴)
- `respiration_C_per_C_h  ∝  C_tot_O2 * r / [ N0 * (exp(r*T_end) - 1) ]`

So the "growth rate" **is** the fit parameter r, and "respiration" is an **explicit
analytic function of r**. Regressing growth on respiration regresses r on a function of r.
R1 was circular by construction, not merely contaminated by a divisor.

There is **no independent growth measurement in this dataset** — no OD or cell-count time
series. Growth is inferred from the curvature of the oxygen trace itself.

## Verdict
R1 (and by extension R2–R6, which all need an independent respiratory constraint) **cannot
be done with these data.** Reviving it requires an independent growth readout (OD or direct
counts over time) measured alongside oxygen, so that growth and respiration are two
measurements rather than two parameters of one curve.

## Correction to earlier advice — this affects Figs 1–2
Earlier in this session I told Ilgaz that Fig 1's decoupling claim was "safe" because the
growth activation energy comes from the rising limb. **That reassurance was premised on
growth and respiration being independent measurements, which they are not.** The
consequences to check:

1. "Growth turns over with temperature while respiration does not" compares two components
   of a single oxygen model, not two independent physiological measurements.
2. The activation-energy gap (E_resp − E_growth) is a contrast between **co-estimated
   parameters** of the same fit.
3. CUE = growth/(growth+respiration) and resp_over_growth are ratios of co-estimated
   quantities; r enters both numerator and denominator.

This does **not** mean Figs 1–2 are wrong. Separating curvature (growth) from slope
(respiration) in an oxygen trace is a legitimate inference, and the two components do carry
different information. But the framing must be "two components of one oxygen model," and
the decoupling needs an explicit check against **fit-induced parameter coupling**.

## Recommended check before submission
- Extract the r–K parameter correlation / joint uncertainty from the per-series fits
  (`fit_coefficients_long.csv`, `fit_metrics.csv` carry the SEs).
- Ask whether the observed growth–respiration relationship (and the E_resp − E_growth gap)
  exceeds what the fit's own parameter coupling would generate under simulation from the
  fitted model with the observed noise level.
- If the gap survives that, Figs 1–2 stand as they are with reworded framing. If it does
  not, the decoupling claim needs to be weakened.

## Status of the respiration programme
Closed. R1 retracted and shown unsalvageable; R2 void (documented in
RESP_R1_R2_CORRECTED.md). No further metabolic-model variants should be built on this
oxygen data.
