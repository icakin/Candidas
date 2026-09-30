# Full re-audit of every quantitative claim (independent re-derivation)

Prompted by unreliable results earlier in the session. Each number below was
**re-derived from the rawest available data with a fresh implementation** and compared
against what was previously delivered. Verdicts: VERIFIED, CORRECTED, or NOT RE-RUN.

## VERIFIED — reproduce exactly, stand as delivered

**Figure 3 (main figure)**
- Growth-rate matrix: reproduces to 5×10⁻⁴ (CSV rounding); NaN/no-growth pattern identical.
  Definition confirmed: `base` = median of per-well growth rate over **28–34 °C**; 36–44 °C
  are single-temperature medians. This matches the caption exactly.
- Isolate counts: growth at 40 °C — *C. auris* 12/12, *C. haemulonii* 1/3,
  *C. duobushaemulonii* 0/2, *C. parapsilosis* 1/3 (pooled non-auris 2/8).
  Detected at 44 °C — 10/12, 1/3, 0/2, 0/3.
- Median highest tested temperature with growth: 44 °C (auris) vs 38 °C (each other species).
- Risk difference **+0.750**, Newcombe 95% CI **[+0.332, +0.929]**, Fisher exact
  **p = 7.22×10⁻⁴**.
- Patristic distances: auris–haemulonii 0.3296, auris–duobushaemulonii 0.3545,
  haemulonii–parapsilosis 1.1663, haemulonii–duobushaemulonii 0.2311.

**Figure 4 (etcGEM counterfactual)**
- Sequence-predicted separations: max |ΔTm| = 0.53 °C, max ΔTopt = 1.00 °C.
- Baseline predicted µ at 44 °C consistent across both result files.
- Required uniform separations, re-derived with an independent script:
  **Tm 32.41 °C** (previously 32.44) and **Topt 15.78 °C** (previously 15.76);
  per-relative values all match; haemulonii and parapsilosis again unreachable via Tm
  within ±35 °C.

**Comparative genomics / HSF1**
- Candidate-gene copy numbers: 5/5 spot-checked families reproduce exactly
  (HSP90, HSP70, catalase, SOD, siderophore transporters).
- HSF1 DNA-binding domain **100 % identical across all four auris clades** under three
  independent domain-window definitions (81, 96 and 106 alignment columns) — robust.

## CORRECTED — errors found, now fixed
1. **Mann–Whitney p-value.** Reported 2×10⁻⁴; correct value is **8.9×10⁻⁴**
   (per-isolate, non-growers scored 0, two-sided). Same direction, still <0.001, but the
   figure was wrong by ~4×. Fixed in FIG3_discord_caption.md and FIG3_caption.md.
2. **Patristic distance range.** Figures/captions said *d* = 0.33–0.36; auris–duobushaemulonii
   is 0.3545, so the correct range is **0.33–0.35**. Fixed in both Fig 3 scripts (Python and R),
   both captions, and the PI email; both figures regenerated.
3. **HSF1 close-relative identity.** Stated 92–95 %; across window definitions it is
   **~90–95 %**. Fixed in CANDIDATE_GENOMICS_SCREEN.md.
4. **R1/R2 respiration analyses — retracted** (documented separately in
   RESP_R1_R2_CORRECTED.md): the respiration predictor was total O₂ divided by a
   growth-derived biomass integral, i.e. the prediction target was inside the predictor.
   With growth-independent oxygen, LOSO MAE 0.1418 → 0.1417 (no gain).

## VERIFIED on a second pass — dynamic sparse search
Re-derived with a fresh implementation (dynamic rerank at every step, beam over pairs and
triples among the top singles, −8 °C, detection 0.05):

| | base µ(40) | 25 enzymes, dynamic | best pair | best triple | reached detection? |
|---|---|---|---|---|---|
| duobushaemulonii Tm | 0.435 | 0.361 | 0.389 | 0.381 | no |
| duobushaemulonii Topt | 0.435 | 0.240 | 0.304 | 0.296 | no |
| haemulonii Tm | 0.588 | 0.499 | 0.550 | 0.550 | no |
| haemulonii Topt | 0.588 | 0.307 | 0.460 | 0.460 | no |

Consistent with the delivered ranges (Tm 0.35–0.60, Topt 0.22–0.40, pairs/triples
0.30–0.55); small upward differences are only because this pass used 25 greedy steps
vs 30, and the trajectory is monotone. Nothing approaches detection. The claim that the
growth optimum is unique is also correct: an LP's optimal objective value is unique even
when multiple flux distributions attain it.

## STILL NOT RE-RUN
- The older **FIG_SUPP_etcGEM** transfer/LOOCV figure numbers (out-of-fold MAE,
  trivial-baseline ratios, qO₂/µ range 3.5–19.5) were not re-audited.
- **Figs 1–2**: growth and respiration are co-parameters of one oxygen fit (see
  RESP_R1_FINAL_VERDICT.md). The decoupling / CUE claims need a fit-coupling check.

## Net effect on the manuscript
Figure 3 and Figure 4 stand, with three numeric corrections (all now applied), none of which
changes any conclusion. The respiration modelling is withdrawn. The two items above remain
unverified and should be re-derived before submission.
