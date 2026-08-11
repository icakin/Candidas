# Claude Code prompt — C2b: write up the N0 findings as a rendered report, and tidy the repository (autonomous)

Run from the project root (`.../Candidas TPC/Candidas`). Two jobs, in order: turn the C1 and C2
findings into a single self-contained technical report rendered to PDF, then consolidate the run
artefacts those two prompts left behind. It runs no analysis, refits nothing, and changes no result.

NOTE TO USER: launch in an auto-approving mode. Quarto and tinytex are already installed (C1). Fast —
the render is seconds and the tidy is file moves. Budget under an hour.

CONTEXT. C1 and C2 have both completed. Between them they produced `reports/REPRODUCTION.md`,
`reports/RUNBOOK.md`, `reports/N0_SENSITIVITY.md`, `reports/figures_N0/` and `reports/tools/`, plus
roughly 1 GB of parallel output trees: `results_C1/` (260 MB) and `results_C2_current/`,
`results_C2_ramp/`, `results_C2_nobp/` (225 MB each), and under
`cauris_etcgem/strains/eci_cauris/outputs/`, `supp_data_C1/` (11 MB), `figure3_data_C1/` and
`_percladefit_C1.npy`. Top level has also accumulated `run_c1.sh`, `run_c2_arms.sh`, `Rplots.pdf` and
`.DS_Store`. Those trees are EVIDENCE for the two reports, not scratch — they must be consolidated
and manifested, not casually deleted.

The headline C2 findings the report must carry:
- `respiration = K/N0` exactly (window length cancels), so `log R = log K - log N_inoc - r*delta`,
  verified to sd 1.3e-14. Median `r*delta` = 0.679 (~2x deflation, max 35.7x), peaking at 32-34 C —
  the growth optimum. It does not add noise; it imprints an inverted growth curve.
- `delta` is instrumental: flat with temperature (58-82 min, R2 = 0.013), unrelated to `r`
  (R2 = 0.002), no taxon effect once glabrata is set aside (p = 0.097). CAVEAT that must appear in
  the report: every set point is above ambient, so |T - T_ambient| is perfectly collinear with T and
  that regression cannot separate them; the transient analysis is what carries the argument.
- Two corrections to the original premise: `tau` does NOT scale with distance from ambient — it is
  17-20 min at every temperature; the 37/51/70 min figures were PEAK TIMES, not time constants. What
  scales is the AMPLITUDE: r(A, |delta solubility|) = +0.996.
- `T_internal [°C]` is recorded in all 72 raw exports (last column, beside `Tm`, `p`, `Salinity`) and
  is discarded by `01_convert_xlsx.R`. It dips 4.76 K below its asymptote at 44 C (t = 10 min),
  2.92 K at 34 C, 0.75 K at 22 C, and is still 0.28 K low at 84 min at 44 C; the 44 C plate asymptotes
  at 43.44 C, never reaching set point. `Tm` is constant at the operator's set point and is what the
  optode's compensation used — so the reported O2 CONCENTRATIONS are compensated at a temperature
  wrong by up to ~4.8 K over exactly the pre-window interval.
- Removing the `r*delta` term destroys the E_R ordering: Clade IV moves from lowest (0.291) to 4th of
  6 (0.597); Spearman between orderings +0.486 (p = 0.33).
- Arms: invariance exactly 0.000e+00 on every growth-side quantity; floor 35/36 tables identical.
  Three of four claims survive. Clade III vs Clade IV becomes 1.00 [0.920, 1.319] under ARM 3, so
  "Clade IV is credibly the cheapest at fever" fails. Fig 4e r(capacity, respiration) goes to
  p = 0.54 — under ARM 3 respiration contains no `r` and capacity is fitted to growth, so that
  correlation exists only via `exp(r*delta)`.
- No endpoint biomass exists anywhere: every column of all 72 raw exports was checked;
  `otu_inoc.csv` is the inoculum, not the end.
- The assay fix is cheap: row D and column 6 are already "No Sensor", i.e. 9 unused sensor positions
  per plate for cell-free blanks, plus one endpoint count per vial.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: reports/
{REPRODUCTION.md,RUNBOOK.md,N0_SENSITIVITY.md}, reports/figures_N0/, reports/tools/,
manuscript/draft/{_quarto.yml,header.tex,nature-communications.csl,references.bib} (to match house
style), scripts/config.R and scripts/01_convert_xlsx.R (for the T_internal claim), and the current
.gitignore. Do NOT re-run any analysis, do NOT refit anything, do NOT change any result.

PART A - the report
- Create `reports/n0/` as a small Quarto project rendering `n0_report.qmd` to BOTH PDF and docx,
  reusing the manuscript's fonts/CSL/header conventions so it looks like it belongs to this project.
  Title it for a reader who has not followed the work: it must stand alone.
- Structure it as an argument, not a log. Suggested sections, adapt as the material demands:
  1. Summary — what was asked, what was found, what should change, in under a page.
  2. The identity — derive `respiration = K/N0` and show the window length cancels; state
     `log R = log K - log N_inoc - r*delta` and the 1.3e-14 verification. Make clear this is a
     structural property of the estimator, not a bug.
  3. What `delta` is — the evidence it is instrumental, INCLUDING the collinearity caveat stated
     plainly rather than buried.
  4. The transient — A, tau and O2_eq against temperature; the amplitude-vs-solubility correlation
     (+0.996); the correction that tau is constant and it is the amplitude that scales.
  5. The recorded internal temperature — the T_internal trace, the numbers above, the compensation
     mismatch, and what it implies for the O2 values themselves. This is a NEW finding: give it its
     own section and be explicit that it is a recomputation available from data already in hand, not
     a new experiment.
  6. The three arms — invariance, the reproducibility floor, and the full comparison table.
  7. What survives — the four claims with SURVIVES / WEAKENS / FAILS and the numbers behind each.
  8. Recommendation — which arm, what to change in the manuscript, and what would settle the
     remaining uncertainty.
  9. What to add to the assay — cell-free wells in the 9 unused positions, endpoint counts, and the
     forward-vs-backward N0 comparison they enable.
  10. Appendix — derivations, the phi/psi identifiability argument, methods, and a table of every
      figure and table the report cites with its path.
- Include the figures from `reports/figures_N0/`. Every number in the text must trace to a file in
  `reports/` or a run tree; add a short provenance line under each table.
- Be even-handed. Three of four claims survived and the headline result strengthened — say so as
  plainly as the failures. Where ARM 3 drives a conclusion, state explicitly that ARM 3 is a known
  lower bound and explain why it is still informative as a DIAGNOSTIC (it isolates the mechanism)
  even though it is not the preferred estimator.
- State clearly which findings are method-level (they apply to `icakin/OxygenModel` and the published
  Cakin et al. 2026 analysis) and which are Candidas-specific.
- Render, verify the PDF opens, has no missing figure and no unresolved cross-reference, and report
  the page count.

PART B - consolidate the run artefacts
- Introduce `runs/` at the project root and move, preserving content:
      results_C1/          -> runs/C1_reproduction/
      results_C2_current/  -> runs/C2_arm_current/
      results_C2_ramp/     -> runs/C2_arm_ramp/
      results_C2_nobp/     -> runs/C2_arm_nobp/
  and under cauris_etcgem/strains/eci_cauris/outputs/:
      supp_data_C1/        -> runs/C1/supp_data/
      figure3_data_C1/     -> runs/C1/figure3_data/
      _percladefit_C1.npy  -> runs/C1/_percladefit.npy
- Update every reference to the old paths — in reports/, scripts/, prompts/, run scripts, .gitignore
  and the two rewritten comment blocks — so nothing dangles. Grep for the old names afterwards and
  show the grep returns nothing but historical prose.
- Write `runs/MANIFEST.md`: for each tree, what produced it (which prompt, which commit, which
  command), when, the env from `env/versions.json`, what it is evidence FOR, and a checksum index of
  its tables. The manifest must make each tree auditable even if it is later pruned.
- PRUNE ONLY WHAT IS PROVABLY REGENERABLE AND UNCITED. Classify every file in the run trees as
  CITED (referenced by a report), EVIDENCE (a table underpinning a claim), or REGENERABLE (figures
  and .rds reproducible by re-running). Keep all tables. Keep anything cited. For the rest, report
  the size saved and delete only after listing what will go. NEVER delete `results/`, `data/`,
  `outputs/supp_data/`, `outputs/supp_data_original_backup/` (that mixture is the evidence for the
  provenance finding), `outputs/figure3_data/`, or anything under `manuscript/`.
- Housekeeping, each verified before acting: move `run_c1.sh` and `run_c2_arms.sh` into `scripts/`
  (or fold them into `run_all.sh` if they are now redundant — check first); delete `Rplots.pdf` and
  any `.DS_Store`; remove any stale LibreOffice `.~lock.*#` files under `manuscript/`; and handle the
  7 stale duplicate manuscript figures C1 found in `results/figures/` — CONFIRM they are byte-identical
  duplicates of the `results/figures/manuscript/` copies before removing, and if they are not
  identical, leave them and report the difference.
- Extend `.gitignore` for the `runs/` convention, keeping tables tracked and excluding bulk
  (.rds, large figure directories). State the policy in one line at the top of that gitignore block.
- Update `README.md` and `reports/RUNBOOK.md` to document the layout: what `results/` is (canonical,
  shipped), what `runs/` is (non-canonical, per-prompt, manifested), and where reports live.

VERIFY (report all)
1. `reports/n0/n0_report.pdf` and `.docx` exist; page count; no missing figure, no unresolved
   reference; every figure in `reports/figures_N0/` either used or explicitly noted as omitted.
2. Every number in the report traces to a named file; list any that does not.
3. Old paths (`results_C1`, `results_C2_*`, `supp_data_C1`, `figure3_data_C1`) return no live
   references anywhere in the repo.
4. `runs/MANIFEST.md` exists and covers every tree with producer, date, env, purpose and checksums.
5. Disk before and after, per tree, with the prune list itemised.
6. `results/`, `data/`, `outputs/supp_data/`, `outputs/supp_data_original_backup/` and `manuscript/`
   byte-identical to their pre-run state — verify by checksum, not assertion.
7. The 7 `results/figures/` duplicates: confirmed identical and removed, or left with the difference
   reported.
8. `bash scripts/run_all.sh --help` (or a dry run) still resolves all paths after the move.
9. README.md and RUNBOOK.md describe the new layout.

CONSTRAINTS
- No analysis. Do not re-run any pipeline stage, refit any model, or alter any number. If writing the
  report reveals a numerical inconsistency between the two source reports, FLAG it in a "discrepancies"
  subsection — do not resolve it by recomputing.
- Move, don't copy-then-delete; nothing may be lost in transit. Every deletion listed before it happens.
- `N0_BACKPROJECT` stays TRUE and `config.R` behaviour stays unchanged.
- The report must be readable by someone who has not followed this work, and honest about what is
  still unresolved — particularly that the T_internal correction has not yet been applied, so the
  fit-window contamination question remains open.
- Autonomous; commit in parts: "C2b: N0 technical report (Quarto -> PDF + docx)",
  "C2b: consolidate run trees under runs/ with manifest",
  "C2b: prune regenerable artefacts, housekeeping, document the layout".
```
