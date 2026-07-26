# Prompts — Candidas project

Each file here is a prompt to paste into **Claude Code** (run from the project root
`.../Candidas TPC/Candidas`). Most are autonomous and modify code, run analyses, and
commit; launch Claude Code in an auto-approving mode (accept edits + allow commands,
e.g. `--dangerously-skip-permissions`) to run them unattended. They build on each
other, so run them in the order below.

Convention follows the `etcGEMs` repo (`.../MICROADAPT/etcGEMs/prompts/`): the `C`
series is the Candidas thermal-performance + etc-GEM application.

Status legend: ✅ run · ▶ in progress · ⏳ pending

## Run order

**Reproducibility**

1. `C1_environment_and_end_to_end_reproduction_prompt.md` — ✅ pinned R/Python/Stan/Quarto, made the
   pipeline headless, re-ran everything into parallel trees and wrote `reports/REPRODUCTION.md`.
   **Outcome:** the pipeline now runs unattended in 49 min. The respirometry half reproduces (worst
   disagreement 2.4×10⁻¹⁰; every Bayesian number within 0.4 % relative). The etc-GEM half does not —
   capacity regenerates 1.33–1.37× higher, traced to an anchor change, and `figure3_data/` uses a
   temperature grid the current code cannot emit. All seven predicted inconsistencies confirmed.

**Respirometry — the load-bearing checks**  *(C3, C4: pending, NOT YET WRITTEN — the
summaries below are a specification, not a prompt you can paste)*

2. `C2_n0_backprojection_decision_prompt.md` — ✅ the N₀ back-calculation, in one prompt.
   Establishes that δ is the vial/optode equilibration time, not a biological lag (δ ≈ 58–82 min at
   every temperature, near-identical across taxa; pre-peak O₂ is a first-order approach to
   equilibrium). Characterises that transient, bounds the two temperature-dependent errors it creates — the thermal-ramp error in
   `exp(r·δ)` and possible contamination of the *fit window* by residual drift, which could inflate
   K at the hot end — then re-runs everything downstream under three arms (current / ramp-aware /
   no-back-projection as a lower bound) and tabulates every load-bearing number side by side. Ends
   by comparing forward against backward back-calculation and specifying the two assay additions
   (cell-free wells at every temperature; one endpoint count per vial) that would settle it without
   disturbing the sealed measurement. Changes no default. *(Merged the former `C2a`, now in
   `archive/`.)*
   **Outcome:** `respiration = K/N₀` exactly, so `log R = log K − log N_inoc − r·δ` (verified to
   sd 1.3×10⁻¹⁴); the term peaks at 32–34 °C and imprints an inverted growth curve rather than adding
   noise. δ confirmed instrumental. Removing the term destroys the E_R ordering (Clade IV: lowest →
   4th of 6). Invariance 0.000e+00; floor 35/36 tables. Three of four claims survive; "Clade IV is
   credibly the cheapest at fever" fails and Fig 4e goes to p = 0.54. **New finding:** `T_internal`
   is recorded in all 72 raw exports and discarded by `01_convert_xlsx.R` — the optode compensated at
   `Tm` while the sample sat up to 4.8 K colder.

3. `C2b_report_and_repo_tidy_prompt.md` — ✅ wrote the C1/C2 findings up as a standalone Quarto
   report rendered to PDF + docx, then consolidated the ~1 GB of parallel run trees under `runs/`
   with a manifest, pruned only what is provably regenerable and uncited, and documented the layout.
   No analysis, no refits, no number changes.
   **Outcome:** `reports/n0/n0_report.pdf` (17 pp) + `.docx`, carrying the φ/ψ identifiability
   argument — φ, the true growth over δ, is *not identified* by the oxygen record, so the three arms
   are three priors rather than three estimates. `runs/` holds the four trees with
   `runs/MANIFEST.md`; 654 MB pruned, every table kept. Of the 7 duplicate figures in
   `results/figures/`, only 5 are byte-identical — two are superseded renders — so all seven were
   left and the difference reported.

4. `handover_and_push_prompt.md` — ✅ wrote `HANDOVER.md` and pushed C1/C2/C2b to GitHub as
   `gyd/c1-c2-review` branches (submodule first), for review rather than for unexamined merge.

5. `C3` — ⏳ the 88 hand exclusions and the hot-end coverage. 43% are flagged by no automatic
   criterion, they are directional (lower `r` in 14/18 testable cells), and *C. parapsilosis*
   has no data above 40 °C yet is tabulated and plotted to 44 °C. Document each exclusion,
   model the non-respiring wells as zeros rather than as missing, and stop extrapolating.
6. `C4` — ⏳ model specification: correlated isolate random effects `(1|i|Isolate)`, a
   replicate/plate level, honest measurement error (`se_r/r` ≈ 0.01 is ~11× smaller than the
   observed replicate reproducibility), Pareto-k on the LOO, and the curved-monotone
   respiration alternative that `08_bayesian_models.R:157-186` says wins but `:187` disables.

**etc-GEM — the load-bearing checks**  *(C5-C8: pending, NOT YET WRITTEN — same)*

7. `C5` — ⏳ build the model PER CLADE. `build()` is currently called once from one MW file
   keyed to the Clade I proteome, so the four a-priori curves are bit-identical by
   construction and "thermotolerance is not sequence-encoded" is untested.
8. `C6` — ⏳ the a-priori curve. It is a flat-topped plateau (the biosynthesis cap), not a
   TPC, and the pool budget is bisected to a measured growth rate — a fourth fitted parameter.
   Make the genome-only prediction a real prediction and report its (negative) R² honestly.
9. `C7` — ⏳ capacity as an amplitude. `kcat_scale` scales growth exactly linearly, so
   corr(kcat_scale, measured peak) = 0.996 by construction and a bare 3-parameter MMRT curve
   reproduces the published R² to 3–4 decimals. Either demonstrate the network contributes, or
   reframe the claim.
10. `C8` — ⏳ Methods provenance: Topt is `0.5·Tm + 10.107`, not TOME; ΔCp is a single global
   constant chosen because the fit demanded it; sector fractions are hard-coded, not from a
   *C. auris* proteome; Tm coverage is 67%, not 47%. Align the text with the files.

## The through-line

Measured O₂ respirometry gives growth and respiration TPCs for four *C. auris* clades and
*C. parapsilosis*; their different thermal sensitivities put the carbon-economy optimum below
37 °C and make a febrile host a quantifiable carbon tax. A sequence-grounded etc-GEM then asks
how much of the clade differences the genome alone explains. The C-series exists to make each
step of that chain reproducible and to test the parts of it that are currently generated by
the setup rather than by the data.
