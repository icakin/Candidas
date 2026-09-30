# Milestone 1 — enzyme pool-binding test (auris ecGEM)

Question (ChatGPT non-negotiable #1): at observed growth rates, does the protein
pool bind (load-bearing) or is it slack (ecGEM collapses to plain GEM)?

Method: auris re-keyed model on YMS (rich scenario, glucose uptake -10, pooled AA);
sMOMENT single-pool constraint  sum_r (MW_r/kcat_r)*(v_fwd+v_rev) <= P.
MW from sequence (enzyme_mw_auris.csv); kcat = EC-class median prior (Bar-Even-style,
1/s: EC1 13.7, EC2 13.7, EC3 79, EC4 15, EC5 6.9, EC6 10); representative-isozyme MW.

Result:
| pool P (g enzyme/gDW) | ecGEM growth /h |
|---|---|
| infinite (plain GEM)  | 2.045 |
| 0.25 (realistic metabolic protein) | 0.758 |
| 0.10 | 0.316 |
| 0.05 | 0.155 |

Measured auris mu_max ~0.80 /h.

Verdict: POOL BINDS (2.7x reduction vs plain GEM at realistic P). Load-bearing.
At a literature protein budget with GENERIC (untuned) EC kcats, ecGEM ~0.76/h ~=
measured 0.80/h -> correct regime. Enzyme layer is viable; proceed to real kcats
(DLKcat) + calibration + the temperature layer.

Caveats: generic EC-class kcats (not species-specific yet - no discrimination until
DLKcat); representative-isozyme MW (OR-GPRs, complexes not summed); single P not yet
calibrated per species. These are refinements, not blockers.
