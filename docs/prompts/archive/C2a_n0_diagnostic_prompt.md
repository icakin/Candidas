# Claude Code prompt — C2a: characterise the equilibration transient and bound the N0 back-calculation error (read-only diagnostic, autonomous)

Run from a SEPARATE WORKTREE of the project (see NOTE TO USER). Read-only diagnostic. It establishes
what `delta` actually is, characterises the physico-chemical equilibration transient that dominates
the pre-window oxygen record, bounds two temperature-dependent errors it creates, and specifies the
measurements that would remove them. It runs NOTHING downstream, refits no Bayesian model, writes to
`reports/` only, and touches no shared state. It is the first half of C2, split out so it can run
while C1 is still going.

NOTE TO USER: launch in an auto-approving mode, in a separate LOCAL CLONE so it cannot collide with
the C1 run's git index. Use a clone, not `git worktree add` — a worktree writes to the main repo's
`.git/worktrees/` and can hit `index.lock` while C1 is committing; a clone only reads:

    git clone . ../candidas-c2a
    cd ../candidas-c2a   # launch here

Everything this prompt needs is tracked (all raw data/*_Oxygen.csv, and results/tables/). The one
exception is `Oxygen_All_Long.csv`, which is gitignored and must be reconstructed — see below.

Cheap: curve fits to short pre-peak segments plus table arithmetic. No MCMC. Safe to run alongside C1,
which will be in its brms stage and wants the cores.

CONTEXT — READ THIS CAREFULLY, IT CORRECTS AN EARLIER MISREADING.

Per-cell respiration is `K / N0`, with `N0 = N_inoc * exp(r*delta)`, so
`log R = log K - log N_inoc - r*delta`.

The cell density in a sealed vial CANNOT be measured at the window start without breaking the seal
and destroying the closed-system oxygen measurement. The back-calculation is therefore a considered
design choice, not an oversight, and the question is which back-calculation is best supported — not
whether to do one.

`delta` is NOT a biological lag. It is the time for the vial and optode to reach physico-chemical
equilibrium after the liquid is added, and the raw traces show this plainly. Pre-peak oxygen follows
a first-order approach to equilibrium whose amplitude and time constant both scale with the distance
from ambient (Clade1 R1: +0.88 mg/L over 37 min at 22 C; +2.78 over 51 min at 34 C; +5.40 over 70 min
at 44 C). And `delta` itself is ~58-82 min at EVERY temperature with no trend, and near-identical
across taxa (Clades I-IV and para all 60-90 min; only glabrata is longer at 100-148 min, consistent
with its known signal-to-noise limit). A biological lag would be far shorter near 34-36 C where
growth is ~3x faster than at 22 C. It is not.

Two consequences follow, and they are the subject of this prompt:

  (i)  THERMAL-RAMP ERROR. During `delta` the sample was climbing from ambient to set point — a ~2 C
       climb at the bottom of the range and ~24 C at the top. The cells were therefore NOT growing at
       the final-temperature rate `r` for the whole interval, so `exp(r*delta)` over-corrects, and it
       over-corrects most at the temperature extremes. This distorts precisely the axis the paper is
       about.

  (ii) FIT-WINDOW CONTAMINATION. If the transient is not fully settled at the oxygen peak, the early
       part of the FIT WINDOW carries residual abiotic drift on top of biological drawdown. That
       would inflate `K`, and inflate it more at high temperature where the transient is larger and
       slower — a candidate mechanism for "respiration rises monotonically to 44 C". There are NO
       cell-free wells anywhere in the design (row D and column 6 are "No Sensor"), so this can be
       bounded from the transient's own shape but not measured directly.

The same construction is in the published method (`icakin/OxygenModel`, `CUE.R:222`), so any finding
here is method-level, not Candidas-specific.

---

```
Work AUTONOMOUSLY end to end; print a summary. Commit ONLY inside this worktree, on a branch named
c2a-n0-diagnostic. Read first: scripts/config.R (the N0_BACKPROJECT block ~446-456, get_baseline(),
resp_model(), the carbon constants), scripts/07_oxygen_fits.R (delta = fit_start_time at ~line 361;
N0 at ~920; respiration and CUE at ~900-1000), scripts/03_trimming.R (how the peak and the window
start are chosen), scripts/04_trim_selector.R:196-260, and
results/tables/{derived_N0_R_results_with_carbon.csv, manual_fit_windows.csv,
Oxygen_Trimmed_Series_Metadata.csv, otu_inoc.csv}.

`Oxygen_All_Long.csv` is NOT committed. Reconstruct it in a scratch directory by reimplementing
02_longdata.R's reshape from data/*_Oxygen.csv, and VERIFY the reconstruction against the derived
table (same (T, OTU, Replicate) keys, same series lengths). Do not write it into results/.

PART A - establish what delta is
- Report `delta` by temperature and by taxon (median, IQR, CV). Test formally whether `delta` depends
  on temperature, on taxon, and on the fitted `r`. A biological lag would shorten sharply near the
  growth optimum; an instrument timescale would not.
- Regress `delta` on |T - T_ambient| (take ambient = 20 C, and state the assumption). Report whether
  the instrumental reading is supported and quantify the residual taxon effect (expect glabrata to
  separate).
- Report the distributions of `r*delta` and `exp(r*delta)` overall, by temperature and by taxon
  (median, IQR, 95th pct, max), and the value at 22, 34 and 44 C per taxon.
- Per taxon compute three Boltzmann-Arrhenius slopes — E_K (raw volumetric K, no N0), E_R as shipped,
  and E_R with the `r*delta` term removed — over the full range AND the rising limb (T <= 34 C).
  Report the artefact (E_K - E_R) and whether the taxon ORDERING of E_R survives removal.
  Frame this as "how much of E_R is the correction", NOT as evidence the correction is wrong.

PART B - characterise the equilibration transient
- For every series, fit a first-order approach to equilibrium to the PRE-PEAK segment:
      O2(t) = O2_eq - A * exp(-t / tau)
  Report per-series A (mg/L), tau (min) and O2_eq, plus fit quality. Handle series that do not rise.
- Show how A and tau vary with temperature and with |T - T_ambient|. Report medians at every set
  point. Test whether tau differs by taxon (it should not, if it is instrumental).
- Report the equilibrium oxygen O2_eq against temperature and check it against the known temperature
  dependence of oxygen solubility in water. Agreement is strong independent evidence the transient
  is physico-chemical; disagreement needs explaining.

PART C - bound the two errors
(C1) FIT-WINDOW CONTAMINATION
- Extrapolate each fitted transient FORWARD past the peak and integrate the residual abiotic drift
  over the actual fit window. Express it as a fraction of the observed biological drawdown over the
  same window, and as an implied fractional bias in `K`.
- Report that fraction by temperature and by taxon. State explicitly whether it grows with
  temperature, and whether it is large enough to account for any part of the monotonic rise in
  respiration to 44 C.
- Be honest about the limits of this: with no cell-free wells the extrapolation is model-based, so
  give a sensitivity band over plausible tau (e.g. the interquartile range of fitted tau at that
  temperature) rather than a single number.

(C2) THERMAL-RAMP ERROR IN exp(r*delta)
- The current correction assumes growth at the final-temperature rate `r` throughout `delta`. Build a
  simple alternative in which the sample temperature follows the fitted equilibration trajectory
  (use tau from PART B as the thermal time constant) and the instantaneous growth rate follows each
  isolate's own fitted thermal performance curve at the instantaneous temperature. Integrate to get a
  ramp-aware N0.
- Report the ratio ramp-aware N0 / current N0 by temperature and taxon, and the resulting change in
  E_R per taxon. State whether the over-correction is largest at the temperature extremes as
  predicted, and whether the taxon ordering of E_R moves.
- This is a BOUND, not a replacement estimator. Present it as "the current correction is wrong by
  about this much, in this direction", not as a new default.

PART D - what would actually settle it
- Set out the two available back-calculations and compare their assumptions explicitly:
    FORWARD (current):  N0 = N_inoc * exp(r*delta)
        assumes exponential growth at rate r through the equilibration/thermal ramp — an interval the
        oxygen record cannot see and during which the cells were not at temperature.
    BACKWARD (proposed): N0 = N_end * exp(-r*T_window)
        assumes exponential growth over the fit window only — which the oxygen model ALREADY assumes.
        Requires one endpoint cell count per vial, taken after the run, so the sealed measurement is
        never disturbed.
- State plainly that the backward route adds no assumption the model does not already make, and that
  running BOTH yields the missing quantity empirically: the ratio
      [N_end * exp(-r*T_window)] / [N_inoc * exp(r*delta)]
  measures how much growth actually occurred during equilibration.
- Check whether any endpoint biomass exists that would allow this NOW — search the repo and
  data/ for endpoint counts, final OD, or flow-cytometry data, and report what is and is not there.
  Note that the OxygenModel paper's flow-cytometry comparison data may permit a retrospective test on
  the 15 bacterial taxa; flag it if so, but do not attempt to obtain it.
- Specify the two additions to the assay, concretely (what to measure, on how many wells, when):
  (1) cell-free medium-only wells at EVERY temperature, to measure the transient directly instead of
      extrapolating it; (2) one endpoint count per vial, to anchor N0 backwards.

PART E - write up
- `reports/N0_DIAGNOSTIC.md` with all four parts, ACTUAL numbers throughout, and figures in
  `reports/figures_N0a/` (delta vs T and vs taxon; A and tau vs T; O2_eq vs T against the solubility
  curve; contamination fraction vs T; ramp-aware/current N0 ratio vs T; E_R per taxon under each
  treatment).
- A RECOMMENDATION section answering four questions directly: (1) is `delta` instrumental — yes or
  no, with the evidence; (2) how much of the reported E_R and its taxon ordering is attributable to
  the correction; (3) does residual equilibration contaminate the fit window, and by enough to matter
  at the hot end; (4) is the two-arm run (C2b) still a decision or now only a robustness check.
- Quote the two contradicting comment blocks (config.R:446-456 and 04_trim_selector.R:247-252)
  verbatim and state which, if either, the evidence supports. Do NOT edit them here — that is C2b.
- State explicitly whether each finding applies to the published method as well as to Candidas.

VERIFY (report all)
1. Reconstructed long table matches the derived table on keys and series lengths.
2. delta vs temperature, taxon and r; the |T - T_ambient| regression; the glabrata residual.
3. r*delta and exp(r*delta) distributions; values at 22, 34, 44 C per taxon.
4. E_K, E_R (shipped), E_R (term removed), full range and rising limb; ordering survives or not.
5. A, tau, O2_eq per temperature; tau by taxon; O2_eq against the solubility curve.
6. Fit-window contamination as a fraction of biological drawdown and as implied bias in K, by
   temperature, WITH a sensitivity band over tau.
7. Ramp-aware / current N0 ratio by temperature and taxon, and the change in E_R per taxon.
8. What endpoint biomass exists in the repo, if any.
9. reports/N0_DIAGNOSTIC.md + reports/figures_N0a/ written; nothing else modified.

CONSTRAINTS
- READ-ONLY with respect to the analysis. Do not modify config.R, any script, results/, results_C1/,
  outputs/, or data/. Write only to reports/ and a scratch directory inside this worktree.
- Do not run 02, 03, 07, 08, 09, 10, 12 or 13. Do not fit any Bayesian model. Do not re-run the
  etc-GEM. Reconstruct what you need and do the arithmetic.
- Do NOT recommend measuring cell density in the vial at the window start. It cannot be done without
  breaking the seal and destroying the measurement. Work within that constraint throughout.
- Treat the back-calculation as a considered design choice. The question is which back-calculation is
  best supported, not whether to back-calculate.
- Stay off the shared git index: commit only on branch c2a-n0-diagnostic inside this worktree.
- Keep it single-threaded / low core count — C1 is running concurrently and needs the CPU.
- Report ACTUAL numbers everywhere, and give ranges rather than point estimates wherever the answer
  rests on extrapolating a transient we have no blanks for.
- Autonomous; commit in parts: "C2a: establish delta as an instrumental timescale",
  "C2a: characterise the equilibration transient (A, tau, O2_eq vs T)",
  "C2a: bound fit-window contamination and thermal-ramp error",
  "C2a: forward vs backward back-calculation + diagnostic report".
```
