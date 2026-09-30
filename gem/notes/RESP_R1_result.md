# ⚠ SUPERSEDED — DO NOT USE

The conclusions in this file are retracted. The respiration predictor contained the
prediction target (O2 divided by a growth-derived biomass integral), so the reported
oxygen signal was leakage. See RESP_R1_R2_CORRECTED.md.

---

# Respiration model — R1 (oxygen-constrained forward): gate result

**Data:** measured per-well specific respiration (`respiration_C_per_C_h`, biomass-normalized
O₂ respiration — an independent measurement, not derived from the growth-rate target) and
growth (`growth_C_per_C_h`), growth-positive wells only, mapped to the four species.

## Two findings

**1. Respiration and growth UNCOUPLE at high temperature.** In every species, measured
O₂ respiration *rises* monotonically with temperature (e.g. C. auris 0.28 → 0.93
C·C⁻¹·h⁻¹ from 22→44 °C) while growth peaks (~32–36 °C) then falls. So cells burn *more*
O₂ per biomass and grow *less* at high T. This overturns the "cells stop respiring at the
boundary" premise: a pure O₂ *upper-bound* constraint cannot suppress predicted growth,
because measured O₂ is high, not low, where growth fails. The informative use of O₂ is as
an energy-spilling / maintenance term (→ R2), not a capacity ceiling.

**2. Raw O₂ carries growth information beyond temperature (out-of-fold).** Leave-one-
species-out, predicting per-well growth rate: temperature-only MAE 0.142 vs temperature +
raw O₂ MAE 0.109 (−23%); O₂ helps for all four species individually (auris 0.151→0.128,
haemulonii 0.175→0.125, duobushaemulonii 0.134→0.095, parapsilosis 0.107→0.088).
Pooled 5-fold: 0.128 → 0.106.

## Implication
The respiration direction is worth pursuing: O₂ is not redundant with temperature. But the
signal is uncoupling (high O₂ + low growth), so the mechanism to add to the etcGEM is a
**temperature-dependent ATP maintenance / energy-spilling term** (R2), calibrated shared
(not per species). **The bar it must beat is the temperature + raw-O₂ regression, LOSO MAE
≈ 0.109.** If the mechanistic etcGEM cannot beat that simple baseline using the same O₂
information, the added mechanism is not earning its place.

## Guardrails honoured
Growth was never used as an input (no O₂/growth, no growth-derived CUE, no
respiration-per-growth). `respiration_C_per_C_h` is normalized by measured biomass (cell
counts), not by the predicted/observed growth rate. Evaluation is out-of-fold (LOSO).

## Next (R2, gated)
Add a shared temperature-dependent non-growth ATP maintenance term to the etcGEM, fit the
shared temperature function on training folds only, predict held-out species/temperatures,
and require it to (i) beat MAE 0.109, (ii) rank species correctly at 40–44 °C, (iii) be
stable to kinetic/threshold choices, (iv) use no growth-derived input. Then R3 (constrained
proteome reallocation), the two-part rate/detection evaluation, and finally inverse accounting.
