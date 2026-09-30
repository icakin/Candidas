# Model-identifiability analysis: is Clade II's deficit enzyme capacity, or a bioenergetic/maintenance cost?

**Goal.** The fitted `kcat_scale` is degenerate — it can absorb metabolic-enzyme abundance, active fraction, saturation, respiratory coupling, non-growth maintenance, or uptake limits. This analysis asks which single (or combined) mechanism, calibrated to **growth only**, best predicts the **held-out** respiration and CUE — i.e. which mechanism the data actually favour for Clade II. It needs no new lab data: your measured per-isolate growth and respiration TPCs plus the existing etc-GEM.

Run in your etc-GEM environment (cobrapy + your `etcgem` provider). Pseudocode below is cobrapy-flavoured; adapt to your model's parameter hooks.

---

## Candidate mechanisms (each gets ONE clade-specific amplitude, exactly like kcat_scale)

| id | mechanism | knob | implementation |
|----|-----------|------|----------------|
| M_E | effective enzyme capacity | `s_E` | global multiplier on the enzyme budget / kcat (your current model) |
| M_M | non-growth ATP maintenance | `s_M` | `NGAM = m0 + s_M·h(T)` — raise the ATPM lower bound |
| M_C | respiratory coupling / P:O | `s_C` | controlled proton leak: add a leak reaction draining the proton-motive force / lowering effective P:O by `s_C`. Do **not** rescale mass-balanced stoichiometry arbitrarily. |
| M_U | carbon-uptake capacity | `s_U` | multiplier on the glucose (carbon) exchange lower bound |
| M_R | respiratory-chain capacity | `s_R` | cap on total ETC / O2 flux `v_ETC_max = s_R·v0` |

**Temperature dependence rule (critical):** each knob is either (a) a *constant clade offset*, or (b) shares one *globally-estimated* temperature function `h(T)` with only the amplitude varying by clade. Never estimate a separate value at every temperature — that guarantees non-identifiability.

---

## The oxygen tie-breaker (do this or respiration predictions are meaningless)

At max growth the O2 flux is generally **not unique**. Use a fixed lexicographic objective for *every* model before reading off respiration:

```python
def predict_growth_and_O2(model):
    sol = model.optimize()                     # 1. max growth
    mu  = sol.objective_value
    with model:                                # 2. fix growth at 0.999*mu
        model.reactions.BIOMASS.lower_bound = 0.999*mu
        # 3. secondary objective: minimise total carbon uptake (or enzyme usage)
        model.objective = model.problem.Objective(
            model.reactions.EX_glc.flux_expression, direction="min")
        s2 = model.optimize()
        O2 = -s2.reactions.EX_o2.flux          # respiration at that vertex
        # 4. report the O2 FVA interval too, so you know if it's well-determined
        lo,hi = flux_variability_analysis(model, ["EX_o2"], fraction_of_optimum=1.0).loc["EX_o2"]
    return mu, O2, (lo, hi)
```
If the O2 interval `(lo,hi)` is wide, a point respiration prediction is not defensible — report the interval and score the observation against it.

---

## Stage A — growth-only calibration, respiration/CUE/fever-tax as pure holdout

For each mechanism M and each clade c:
1. **Fit** the single knob `s` to that clade's measured growth TPC only (same procedure you used for kcat_scale; differential evolution is fine).
2. **Predict** the respiration TPC (via the lexicographic routine above), then CUE = μ/(μ+R) and fever tax, using the growth-fitted `s`. These are never used in fitting.
3. **Score** each mechanism by held-out fit to respiration and CUE (RMSE / R² on log-respiration; and the fever-tax ordering). The winner is NOT the best growth R² (all five will fit growth) — it is the best *held-out respiration/CUE*.

**Discriminating predictions (what each mechanism should do):**
- M_M (maintenance) / M_C (coupling): lower growth **and disproportionately higher respiration per unit growth**, strongest at high T. ← matches your observed growth-down/respiration-up pattern.
- M_U (uptake): lower growth but struggles to raise respiration per growth.
- M_R (ETC cap): lowers respiration and growth together (lower absolute O2) unless a second inefficiency is added.
- M_E (enzyme): reproduces growth amplitude; respiration prediction weak/ambiguous.

So the a-priori expectation is that **M_M or M_C beats M_E on held-out respiration** — that is the whole test.

---

## Stage B — joint growth+respiration fit + hybrids

Only after Stage A, fit the plausible models jointly:
```
log μ_ict ~ N(log μ_pred, σ_μ)
log R_ict ~ N(log R_pred, σ_R)
```
Use the O2-model posterior uncertainty on μ and R where you have it (don't treat extracted means as exact). **Do not** add CUE and fever tax as extra likelihood terms — they are algebraic functions of μ and R; that double-counts. Use them only as posterior-predictive checks.

Then compare a *small* set of hybrids with strong regularisation:
```
M_E+C, M_E+M, M_C+M     with  log s_E, log s_C ~ N(0, 0.25^2)
```
Key readouts:
- coupling alone explains Clade II, and `s_E → 1` once coupling is fitted → **it was bioenergetic, not enzyme amount**;
- or `s_E` stays ~0.4 and `s_C` ~baseline → **enzyme/allocation after all**;
- or a **ridge** (low `s_E` ⇄ poor `s_C` interchangeable) → *the data cannot identify the mechanism* — itself a clean, publishable result and the direct justification for absolute proteomics.

---

## Cross-validation & the real test

- **CV unit:** leave-one-**isolate**-out (or leave-one-temperature-block-out). Never treat replicate–temperature points as independent — errors are nested and curve params couple temperatures.
- **The decisive test:** fit growth on 11 isolates, predict the **full respiration curve of the 12th**.
- **Within-Clade-II ordering:** report whether each mechanism reproduces the 0.37 < 0.47 < 0.51 capacity gradient across the three Clade II isolates. A 4-clade correlation is only n_eff = 4 (smallest permutation p ≈ 0.04); the within-clade gradient is where the real discriminating power is.

## Decision summary
- M_E wins held-out respiration → report capacity as (proxy for) metabolic-enzyme deployment.
- M_M/M_C wins → reframe the clade difference as a maintenance/coupling (bioenergetic) trait; capacity was laundering it.
- Ridge → state that growth+respiration alone cannot identify it; absolute spike-in proteomics (+ cell size, ATP, O2 in the same cultures) is required to partition the branches.
