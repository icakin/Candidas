# Figure 4 — locked state

Frozen 2026-09-04 at commit `ce01596`. Every number the figure prints is listed below with
the script that produces it, so a later edit that changes one of these is visible as a
disagreement with this file rather than as silent drift.

Regenerate the whole figure with:

    python3 gem/17_build_measured_tpc.py     # panel A's measured curve
    python3 gem/26_fig4_tables.py            # flatten the etcGEM JSON into gem/tables/fig4/*.csv
    Rscript scripts/19_fig4.R                # the figure (R, base graphics, since 2026-09-06)

## What the figure asserts

| quantity | value | produced by |
|---|---|---|
| required *T*~m~ separation, fitted formulation | 32.54 °C | `gem/19_etcgem_counterfactual.py` |
| required *T*~m~ separation, corrected unfolding width (w = 2.2) | 18.0 °C | `gem/22_thermal_sensitivity.py` |
| required *T*~opt~ separation | 15.62 °C | `gem/19_etcgem_counterfactual.py` |
| paired Δ*T*~m~, *C. auris* − *C. haemulonii* | 0.411 °C [0.19, 0.63], n = 432 pairs, p = 2.2 × 10⁻⁴ | `gem/24_paired_dedup_audit.py` |
| paired Δ*T*~m~, − *C. duobushaemulonii* | 0.326 °C [0.08, 0.57], n = 440, p = 0.011 | same |
| paired Δ*T*~m~, − *C. parapsilosis* | −0.207 °C, n = 389, p = 0.20 | same |
| paired Δ*T*~opt~, largest | 0.548 °C [−0.00, 1.10], *C. haemulonii*, p = 0.051 | same |
| fold gap, *T*~m~ | 79× (fitted), 44× (corrected width), 20× (measured benchmark) | derived |
| fold gap, *T*~opt~ | 29× | derived |
| measured congeneric benchmark | 1.6 °C over 827 pairs (parental context) | Walunjkar et al. 2025 |
| model-predicted thermal limit | 52.7–54.5 °C across the four species | `gem/22_thermal_sensitivity.py` |
| held-out prediction at each relative's observed failure temperature | 0.79–0.85 of the model's own peak | `gem/25_loocv_nested.py` |

## Three decisions that are easy to undo by accident

**Pairing is by unique protein pair, not by reaction.** Merging the kcat tables on
`reaction` gives 1041 rows from 432 unique pairs, one gene appearing under 23 reactions.
That over-states n AND weights the mean by promiscuity, since highly connected enzymes are
central-metabolism enzymes rather than a random sample. `_paired` in the figure script
collapses to unique (auris gene, relative gene) pairs. Reverting this returns 0.518 °C and
p = 1 × 10⁻¹¹, and makes *C. duobushaemulonii*'s *T*~opt~ difference look significant when
it is not (p = 0.008 reaction-level against p = 0.57 deduplicated).

**The measured TPC counts dead wells as zeros, with the denominator from the raw traces.**
Building it from `derived_N0_R_results_with_carbon.csv` instead drops temperatures where
every well died, because a flat trace cannot be fitted, and the species curve then stops
exactly where the collapse it is meant to show occurs.

**The red cross in panel A means every well died, not a zero median.** Five of fifteen
*C. haemulonii* wells still grow at 40–44 °C. Only *C. duobushaemulonii* (40–44) and
*C. parapsilosis* (42–44) earn a cross.

## Known limitation, not yet resolved

Seq2Tm feeds padded batches to the model with no attention mask, so a protein's predicted
*T*~m~ depends on which proteins share its batch: batch=4 against batch=1 moves the same
protein by up to 5.36 °C, sd 1.27 °C, on 200 *C. auris* enzymes. `gem/tables/thermal_tm.csv` was
produced at batch=4 and is reproduced exactly at that setting (max difference 0.0000 °C,
r = 1.000000, verified in file order). The paired means should be robust, since 432 pairs
average unbiased noise and paired proteins sit in different files, but this has not been
confirmed: `gem/15_run_seq2tm.py --seqs-from ... --batch-size 1` followed by
`gem/audits/padding_sensitivity.py` tests whether the artefact is biased per species. Until that
runs, the methods should state the limitation.
