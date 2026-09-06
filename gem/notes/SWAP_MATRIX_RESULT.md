# kcat swap-matrix result (auris network, scale-matched profiles)

Method: DLKcat kcats per reaction per species; each profile normalized to a common
geomean (10/s) to isolate PATTERN from SCALE (addresses the kcat/P confound). Enzyme
demand = sum_r meanMW_r * |v_r| / kcat_r for the auris pFBA/FBA flux state, under each
profile. Permutation = auris kcat values shuffled across reactions.

| profile             | enzyme demand (g/gDW) | rel. to auris own |
|---------------------|----------------------:|------------------:|
| auris (own)         | 6.07 | 1.00 |
| haemulonii          | 6.10 | 1.005 |
| duobushaemulonii    | 5.74 | 0.95 |
| parapsilosis        | 3.52 | 0.58 |
| PERMUTED (auris)    | 29.8 | 4.91 |
| flat (no pattern)   | 1.79 | 0.29 |

## Interpretation
1. SIGNAL, not noise: permuting the kcat pattern raises demand ~5x. The kcat pattern
   is genuinely aligned with metabolic flux demand -> the 6.7x cross-species spread is
   not pure DLKcat noise.
2. Signal is CONSERVED among close relatives: auris ~= haemulonii ~= duobushaemulonii
   (within 5%); only distant parapsilosis differs (0.58x). => baseline kcat will NOT
   strongly discriminate the held-out close relatives (hae/duo) from auris. Species
   discrimination of hae/duo TPCs must come mainly from the THERMAL layer (Topt/Tm),
   consistent with the review's expectation.

## Caveats
- Single FBA flux state (alternate optima exist); permutation contrast is robust to
  this but the auris~hae closeness is evaluated on auris's own flux state.
- DLKcat ~1-log error remains; permutation test argues the pattern exceeds noise, but
  fine per-reaction values are still uncertain -> use as prior, ensemble later.
