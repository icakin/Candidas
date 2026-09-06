# ⚠ SUPERSEDED — DO NOT USE

The conclusions in this file are retracted. The respiration predictor contained the
prediction target (O2 divided by a growth-derived biomass integral), so the reported
oxygen signal was leakage. See RESP_R1_R2_CORRECTED.md.

---

# Respiration model — R2 (oxygen + temperature-dependent maintenance): result

**Design (as specified):** constrain to measured O₂; add a non-growth ATP maintenance
term; infer a SHARED temperature-dependent maintenance function on training species only;
freeze; predict fully held-out species (LOSO). Maintenance parametrized as rising from a
threshold, m(T)=m₁·max(0,T−T₀), and (separately) exponential; the etcGEM's O₂→ATP→growth
coupling is the energetic backbone (qO₂ = a·growth + m(T) + spill).

## Result: R2 FAILS the gate — and fails structurally
LOSO per-well growth MAE:
- temperature only: 0.142
- temperature + raw O₂ (regression, the bar): 0.109
- **maintenance energetics: 0.329** — worse than temperature alone.

It is not a tuning failure. The maintenance-energetics form makes predicted growth
**increase** with respiration (respiration is the ATP source that funds growth), but the
data show respiration is **highest exactly where growth is lowest** (Fig R2A: the same O₂
value maps to high growth at low T and low growth at high T — a temperature hysteresis).
Maintenance m(T) can only lower the intercept; it cannot invert the sign of the O₂→growth
relationship. So the model systematically **over-predicts** growth at high temperature
(Fig R2B: high-T points fall well above the identity line), and no shared maintenance
function fixes it.

## Interpretation (this is the informative part)
Oxygen carries genuine phenotype information (R1: it improves out-of-fold growth
prediction as a statistical covariate, 0.142→0.109). But that information is **not
metabolic-energetic** and the etcGEM cannot translate it: within a growth-maximizing
metabolic model, all O₂-derived ATP is available for growth, so rising respiration can
only mean more growth or more maintenance — never the observed "more respiration, less
growth." The rising O₂ at high temperature reflects **uncoupling / futile or dissipative
processes / damage-repair**, which have no representation in flux-balance metabolism.

## Consequence for the remaining planned models
R3 (constrained proteome reallocation) and the inverse-accounting model operate inside the
same growth-maximizing metabolic framework, so they inherit the same structural limit:
they can restrict or re-price capacity, but none represents respiration that rises while
growth falls. They are therefore unlikely to convert the etcGEM into a growth predictor.
The honest reframing is the **two-part outcome** ChatGPT flagged: the etcGEM addresses
metabolic *capacity/rate*, not *establishment/viability*. The missing layer is viability
and stress/damage physiology, not metabolic capacity, maintenance energetics, or proteome
allocation.

## Bottom line
- R1: O₂ has out-of-fold growth information (green light to test mechanism).
- R2: the metabolic model cannot translate it — respiration and growth uncouple in the
  direction energetics cannot represent. **Informative negative.**
- Combined with Fig 4 (enzyme thermostability insufficient), this triangulates the same
  conclusion from a second, independent direction: the determinant of C. auris
  thermotolerance lies outside metabolic-network capacity — in regulation, viability, and
  stress/damage physiology.

## Recommendation
Do not keep building metabolic-framework variants (R3, inverse) expecting a predictor —
the structural result says they will not translate the O₂ signal either. Fold R1+R2 into
the paper as a compact second-line supplement reinforcing Fig 4's conclusion, and (if
mechanism is wanted) pursue the wet-lab directions (thermal proteome profiling / stress
proteomics), which measure the missing layer directly.
