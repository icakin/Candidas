# UNBLINDING PLAN — hashed BEFORE blind_key.csv is opened
# Inputs frozen: blind_states.csv 8ece579650a9458c1b65969642fd3c7429ba9aea02ba2afee6509e8bf810dba9
#                classifier.py    7baa22e6e6cf498abf00b094801b67d4f2c9f9fb57affb4cee2a17663af1547e
#                PROTOCOL_v4.md   84b741f3a2c5ef66c60b39201bba8c58651c6ef01c52cdf286d7ccc77debb4c7

## CLASSIFIER STATUS — RECORDED AS QUALIFIED, NOT AS PASSING

v4 is accepted as a QUALIFIED classifier. It did NOT pass all acceptance criteria.

  * A5 FAILED as originally specified: pooled simulated ambiguity 33.3% against the
    25% criterion.
  * The failure is concentrated in the deliberately dense borderline region
    rt in (0, 0.5), where AMBIGUOUS is the appropriate conservative outcome
    (84.4% ambiguous there).
  * Among non-borderline simulations (r = 0 or rt >= 1), ambiguity was 4.9%.
  * A1, A2, A3, A4, A6 passed: 0.000, 0.983, 0.000, 0.000, 0.000.
  * NO threshold, state definition or real-well classification was changed in
    response to the A5 result.
  * A5 was NOT redefined, reweighted or replaced post hoc. No v5 was created to
    make the paperwork pass. A reweighted simulation may later be reported as a
    descriptive sensitivity analysis; it will not replace this A5 result.

## OUTPUTS TO BE PRODUCED, ONCE, WITHOUT MODIFYING ANY CLASSIFICATION

  T1  state counts by taxon x temperature (all 1,200 non-glabrata wells)
  T2  state counts by isolate x temperature
  T3  proportions with uncertainty clustered at ISOLATE level (cluster bootstrap,
      2000 resamples of isolates within taxon)
  T4  complete 37-44 C table (this dataset has 38, 40, 42, 44)
  T5  confusion matrix: automatic state vs legacy manual include/exclude, and vs
      the legacy growing/respiring/inert labels
  T6  location of all 161 ambiguous wells by taxon, isolate and temperature
  T7  ISOLATE-level transition summary (an isolate is scored by the modal state of
      its 5 wells; ties -> ambiguous)
  T8  sensitivity bounds: ambiguous counted (i) as growth, (ii) as respiration-only
  T9  frozen results file reprocess/unblinded_results.csv, hashed, written BEFORE
      any interpretation

## THE ONE QUESTION TO BE ANSWERED FIRST

Does the Figure 5 pattern - relatives entering respiration-without-detectable-growth
at 40 C while C. auris does not - survive under:
   (1) the primary four-state classification;
   (2) all-ambiguous-as-growing;
   (3) all-ambiguous-as-respiration-only?

Decision rule, fixed now: if the C. auris-versus-relatives contrast holds at the
ISOLATE level under BOTH ambiguity bounds, the result is robust. If it holds under
only one allocation, it is reported as UNRESOLVED. No other outcome is available.

## NOT DONE AT THIS STAGE
Figures 1-3 are not refitted. No biological contrast beyond the above is
interpreted. The legacy pipeline, removed_points.csv and all legacy outputs remain
untouched.
