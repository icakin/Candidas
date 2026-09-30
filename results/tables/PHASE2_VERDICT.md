# PHASE 2 VERDICT (CORRECTED, FROZEN) — state-transition framework
# Frozen v4 states (classifier sha256 7baa22e6…), never recomputed. Ambiguity handled only
# as two bounds. Isolate is the unit of inference. No window or threshold was optimised.
# Supersedes the first draft, which imputed censored values and over-generalised criterion 2.

## CRITERION 1 — sharp increase in respiration without confirmed growth near 40 °C: MET

P(respiration without confirmed growth), change 38 -> 40 °C:
  relatives median +0.40 (5/8 isolates jump >= 0.4);  C. auris median +0.00 (1/12).
Seven of eight relative isolates reach P = 1.00 at 40 °C; Hae_1724 is the exception (0.00).

### Confirmed-growth loss, CENSORING-AWARE (corrected)
Right-censored isolates are reported as >44 °C and are NEVER imputed as 46 °C.
Comparison by log-rank statistic with EXACT permutation (20,000 permutations, n = 20).

  bound R (ambiguous = respiration-only)
     all isolates          chi2 = 5.19   p = 0.0242    [C. auris 12 (1 censored), rel 8 (1 censored)]
     excluding Hae_1724    chi2 = 15.57  p = 0.0002
  bound P (ambiguous = productive)
     all isolates          chi2 = 3.37   p = 0.0731    NOT significant at 0.05
     excluding Hae_1724    chi2 = 5.42   p = 0.0221

**CORRECTION TO THE FIRST DRAFT.** Imputing censored isolates at 46 °C inflated these
tests (previously reported p = 0.0065 and p = 0.0330). With censoring handled properly the
contrast is significant under bound R, but only MARGINAL under bound P when all isolates
are included (p = 0.073). The claim "significant under both bounds" is withdrawn; it holds
under both bounds only when Hae_1724 is excluded.

Observed losses, no imputation:
  bound R  C. auris  38, 42, 42, 42, 42, 42, 44, 44, 44, 44, 44 °C  + 1 isolate >44
           relatives 36, 38, 38, 38, 40, 40, 40 °C                  + 1 isolate >44
  bound P  C. auris  42, 42, 44, 44, 44 °C                          + 7 isolates >44
           relatives 38, 40, 40, 42, 44, 44 °C                      + 2 isolates >44

### Respiration loss
Detectable respiration never fell below 50% in ANY isolate at ANY temperature; only 1 of
1,200 wells is "no detectable respiration". T_resploss > 44 °C for all 20 isolates.
The respiration-without-confirmed-growth interval is therefore LOWER-BOUNDED
(>= 44 - T_growthloss), never estimated, and its width is not recoverable from this assay.

## CRITERION 2 — substantial oxygen consumption (REWRITTEN): PARTIALLY MET, taxon-specific

Fixed prespecified window 90-390 min (the manuscript's 1.5 h equilibration exclusion plus
EARLY_H = 5 h from scripts/20_fig4.py), identical for every well. Oxygen-floor censoring:
confirmed-growth wells 45.4%, respiration-only 7.0%.

  * **C. parapsilosis retains substantial oxygen consumption after confirmed growth is
    lost** (3.38 mg/L median at 40 °C, against a 1.0 mg/L respiration gate).
  * **C. haemulonii (1.04) and C. duobushaemulonii (1.23) respiration-only wells sit only
    slightly above the detection gate** and cannot be described as substantial.
  * **Across relatives, respiration is reduced to approximately 41-75% of the same
    isolate's confirmed-growing value at 38 °C** (Hae_1768 0.41, para_2051 0.43,
    para_2053 0.46, Duo_1771 0.45, para_2052 0.75).

DIRECTION OF THE CENSORING BIAS: growing wells reach the oxygen floor inside the fixed
window far more often than respiration-only wells, so their consumption is UNDERSTATED.
The true respiration-only / growing ratio may therefore be LOWER than observed, not higher.

WITHDRAWN: "undiminished", "highest rate in the series", and any universal claim that
respiration remains quantitatively unchanged after growth is lost. None is supported.

## CRITERION 3 — higher productive-respiration fraction under both bounds: MET

DEFINITION.
  Numerator   : O2 consumed (fixed 90-390 min window) summed over wells of ONE isolate at
                ONE temperature whose frozen state is 'confirmed growth'
                (upper bound: 'confirmed growth' OR 'ambiguous').
  Denominator : O2 consumed summed over ALL wells of that same isolate and temperature,
                whatever their state.
  Aggregation : per isolate x temperature, then summarised ACROSS ISOLATES by median and
                FULL RANGE. Pooled oxygen-weighted values are secondary description only.

### At 40 °C, per isolate (the median alone is NOT sufficient)
  C. auris   lower bound median 1.00, RANGE 0.15-1.00 ; upper bound median 1.00, range 0.15-1.00
  relatives  lower bound median 0.00, RANGE 0.00-1.00 ; upper bound median 0.59, range 0.21-1.00

RECONCILIATION of the C. auris median of 1.00 with its 5 respiration-only and 3 ambiguous
wells at 40 °C: those eight wells are NOT spread across isolates. They sit in four:
  Clade2_2073  0.15 / 0.15  (4 respiration-only, 1 growth)
  Clade1_2069  0.62 / 0.79  (1 respiration-only, 1 ambiguous, 3 growth)
  Clade2_2072  0.76 / 1.00  (1 ambiguous, 4 growth)
  Clade1_2070  0.81 / 1.00  (1 ambiguous, 4 growth)
The remaining eight C. auris isolates have 5/5 confirmed growth. The median of 1.00 is
driven by those eight and CONCEALS Clade2_2073, which behaves like a relative.

### Pooled oxygen-weighted fraction at 40 °C (secondary)
  C. auris                lower 0.87  upper 0.92   (462 mg/L over 60 wells)
  relatives               lower 0.33  upper 0.67   (120 mg/L over 40 wells)
  relatives, no Hae_1724  lower 0.00  upper 0.50   ( 81 mg/L over 35 wells)

### Collapse 38 -> 40 °C
  lower bound: C. auris median +0.00 (0/12 drop > 0.5); relatives +0.41 (3/8)
  upper bound: C. auris median +0.00 (1/12);            relatives +0.40 (2/8)
  excluding Hae_1724: relatives unchanged (+0.41 / +0.40)
Because growing wells are censored more often, the productive fraction is itself a lower
bound wherever censoring occurs; the C. auris advantage is at least as large as shown.

## STATUS
Criteria 1 and 3 met. Criterion 2 met only for C. parapsilosis; marginal for C. haemulonii
and C. duobushaemulonii. Criterion 1 is significant under bound R but marginal under bound
P with all isolates included.
Outcome label (A or B) is NOT assigned until Phase 3 is complete.
