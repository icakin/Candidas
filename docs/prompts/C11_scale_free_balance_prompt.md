# Claude Code prompt — C11: establish the CUE-optimum claim without the absolute scale (autonomous)

Run from the project root (`.../Candidas TPC/Candidas`). Analysis on existing outputs. It fits no new
oxygen curves, changes no estimator, alters no committed number and edits no manuscript text. It asks
whether the manuscript's headline claim can be stated without depending on N0, the cell-carbon quota,
the respiratory quotient or any conversion constant — and if so, by how much the model-free estimate
differs from the fitted one.

NOTE TO USER: launch in an auto-approving mode. Branch off `gyd/figure-sourcing` (the top of the
current Candidas stack) and open a PR against it. Cheap — regressions and arithmetic on tables that
already exist, plus a bootstrap.

CONTEXT — WHY THIS MATTERS FOR THE MANUSCRIPT.

The headline claim is that the carbon-use-efficiency optimum lies below human body temperature in
every taxon (31.3-34.1 C, P >= 0.999). It is currently established through a fitted CUE curve, which
inherits the N0 back-projection, the per-cell carbon quota, the O2-to-C conversion and RQ = 1 — all
of which are unresolved to varying degrees.

Work in the OxygenModel repository (report D8) established that the claim does not need any of them.
With `G = r * q` and `R = c * K`:

    CUE = G / (G + R) = 1 / (1 + (c/q) * K/r)

Maximising CUE means MINIMISING `K/r`, and the constant `c/q` drops out of the argmin entirely. So
`T_opt(CUE)` is identifiable from the fitted parameters alone. Only the LEVEL of CUE depends on the
conversions. If that holds here, the manuscript can assert its headline UNCONDITIONALLY, and the
whole absolute-scale discussion becomes a footnote about CUE's magnitude rather than a threat to the
main result.

TWO DIFFERENCES FROM D8 THAT MUST BE CHECKED, NOT ASSUMED.

- D8 found the quantity is `(K * O2ref)/r`, not `K/r`, because the OxygenModel pipeline fits a
  NORMALISED trace and `O2ref` is solubility, which falls ~1.36x across its temperature range. The
  Candida pipeline sets `FIT_TO_SPLINE <- FALSE` and fits RAW `Oxygen`, so `K` should already be
  volumetric in mg/L/min and no `O2ref` term should be needed. VERIFY THIS from `07_oxygen_fits.R`
  and `config.R` before relying on it. If K is not what you expect, the balance quantity changes and
  everything downstream with it.
- D8 found a 7.09 C gap between the model-free optimum and the pipeline's fitted one, of which ~3.2 C
  came from the Sharpe-Schoolfield functional form alone and 2.51 C from a delay override. Candida
  reports 31.3-34.1 C against a 37 C threshold, so a few degrees matters in either direction. A
  downward shift strengthens the claim; an upward one narrows the margin.

The mass-balance half of D8 does NOT transfer: Candida has inoculation density but no endpoint
biomass, so the independent CUE check is not available. State that explicitly so nobody goes looking.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: scripts/config.R (the
carbon constants, N0_BACKPROJECT, FIT_TO_SPLINE, resp_model); scripts/07_oxygen_fits.R (how r, K, G,
R and CUE are formed, and in what units); scripts/09_bayesian_models.R and 12_carbon_tax.R (how
T_opt(CUE) is currently obtained); scripts/08_temperature_equilibration_sensitivity.R and
18_n0_treatment_panel.R (Ilgaz's three N0 treatments — reuse them, do not reimplement);
results/tables/{derived_N0_R_results_with_carbon.csv, fit_coefficients_wide.csv, fig_values.csv,
cue_quadratic_fit_coefs.csv}; and reports/N0_SENSITIVITY.md. Branch:
`git switch -c gyd/scale-free origin/gyd/figure-sourcing`.

PART A - establish what K is, before anything else
- Confirm from the code that the model is fitted to RAW oxygen and that `K` is therefore volumetric
  (mg O2 / L / min) rather than a normalised rate. Report the units explicitly.
- If K IS normalised, work out the correct balance quantity as D8 did and use that throughout,
  reporting the change. Do not proceed on the assumption.
- Derive and state the balance quantity for this pipeline. Show the algebra, including which
  constants cancel in the argmin and which do not.

PART B - the scale-free quantities, per taxon
Work from the fitted r and K only. Do not convert to per-cell carbon anywhere in this part.
- Per isolate and per taxon, fit Boltzmann-Arrhenius to `r(T)` and to `K(T)`, and a unimodal
  alternative where the data support it. Report `E_r`, `E_K` and `E_K - E_r` with intervals.
- Compute the balance quantity against temperature and locate its minimum — the model-free
  `T_opt(CUE)`. Report per taxon with a bootstrap interval over replicates.
- Report the fold-change in the balance across the measured range.
- THE HEADLINE TEST: for each of the five taxa, is the model-free optimum below 37 C, and by how
  much? Report the margin and the bootstrap probability. This is the claim, established without any
  conversion constant.
- State explicitly which quantities carry NO dependence on N0, the carbon quota, the O2-to-C factor
  or RQ. That statement is the point of the prompt.

PART C - compare against the published estimate, and decompose any gap
- The manuscript reports `T_opt(CUE)` of 31.3-34.1 C per taxon. Compare the model-free estimate
  against it, per taxon.
- Where they differ, decompose the gap into its sources, as D8 did:
  * the N0 back-projection (use Ilgaz's existing treatments in 08/18 rather than writing new ones);
  * the functional form — the fitted CUE curve versus a model-free argmin;
  * any solubility or unit term, if PART A finds one applies.
  Report each contribution in degrees, and how much of the total is the fitted curve rather than the
  data.
- If the model-free optimum is LOWER, the claim strengthens and the margin to 37 C widens; say so.
  If it is HIGHER, report the narrowed margin plainly and state whether the claim still holds for
  every taxon. Do not soften either outcome.

PART D - what the manuscript can say
- One paragraph, written to be quoted: which of the manuscript's claims are secure without resolving
  the absolute scale, and which are not. At minimum address: the sub-37 C CUE optimum; "growth is
  more temperature-sensitive than respiration"; the between-taxon ordering; and the absolute CUE and
  per-cell respiration values.
- A short, concrete recommendation on reporting — which quantities should lead, and which should
  carry a stated uncertainty rather than a point value. Do NOT edit the manuscript; this is a
  recommendation for the authors.
- Note explicitly that the mass-balance cross-check available in OxygenModel is not available here,
  and what would make it available (endpoint biomass per vial).

PART E - report
- `reports/C11_scale_free/C11_scale_free.qmd` rendered to PDF, self-contained, with a provenance
  footer mapping every number to its file, figures in its own `figures/`, and every input it reads
  tracked so it rebuilds from a clean clone. If the repository has no report template, follow the
  layout of `reports/n0/` and say so.
- Figures: `E_r` and `E_K` per taxon; the balance quantity against temperature per taxon with its
  minimum marked and 37 C indicated; model-free versus fitted `T_opt(CUE)` per taxon; the gap
  decomposition.
- Be even-handed. If the model-free estimate materially weakens the claim for any taxon, that is the
  finding and it leads.

PART F - PR
- Push and open a PR against `origin/gyd/figure-sourcing` (retarget to `main` once the stack merges).
  Body: the balance quantity for this pipeline, the per-taxon model-free optima, the comparison with
  the published values, the gap decomposition, and an explicit statement that no fitted number in the
  pipeline changes. Do not push to `main` or to any branch in the stack.

VERIFY (report all)
1. The units of K, and the balance quantity derived for this pipeline, with the algebra.
2. `E_r`, `E_K`, `E_K - E_r` per taxon with intervals.
3. Model-free `T_opt(CUE)` per taxon with bootstrap intervals; the margin to 37 C; the probability
   below 37 C for each of the five taxa.
4. The explicit list of quantities free of N0, the quota, the O2-to-C factor and RQ.
5. Model-free versus published `T_opt(CUE)` per taxon, and the gap decomposed by source in degrees.
6. Fold-change in the balance across the range.
7. The one-paragraph verdict on which claims are secure.
8. Report renders from a clean clone; page count; inputs tracked.
9. `results/`, `data/`, `config.R` and the submodule pointer byte-identical, by checksum; no
   manuscript file touched.

CONSTRAINTS
- No estimator change, no constant change, no committed number altered, no manuscript edit.
- PART B must not convert to per-cell carbon anywhere. If a quantity needs the quota or N0 it belongs
  in PART C, and say so.
- Reuse Ilgaz's N0 treatments from `08` and `18`. Do not reimplement them — a second implementation
  is a second thing to keep correct, and D5-D7 in the sister repository showed how easily an
  apparent finding turns out to be an artefact of a reimplementation.
- Verify PART A before relying on it. The equivalent assumption was wrong in the sister repository.
- Report ACTUAL numbers everywhere, with intervals wherever a bootstrap is cheap.
- Autonomous; commit in parts: "C11: establish the balance quantity for this pipeline",
  "C11: scale-free E_r, E_K and the model-free CUE optimum per taxon",
  "C11: compare against the published optimum and decompose the gap",
  "C11: report (Quarto -> PDF) + PR".
```
