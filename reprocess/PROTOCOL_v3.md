# PROTOCOL v3 — written after v2 FAILED simulation calibration
# Labels remain sealed. blind_key.csv has still NOT been read.
# Parents: v1 f9aeed8f...  v1.1 990320e3...  v2 (see PROTOCOL_v2.sha256)

## V2 FAILURE RECORD (sim_results_v2.csv, 1800 traces)
  A1 0.000 <= 0.05 PASS | A2 0.790 >= 0.80 FAIL | A3 0.000 <= 0.50 PASS
  A4 0.000 <= 0.10 PASS | A5 0.422 <= 0.25 FAIL | A6 0.000 <= 0.05 PASS

## DIAGNOSIS

v2 fixed the admission gate and produced a sharp, correct detection boundary at the
intended place (P(growth) = 0.000 below rt 0.5, 0.90 immediately above). Two residual
problems, of DIFFERENT kinds:

  DEFECT C (real, repaired here). At rt_true >= 1 the sole failure mode is reason
  PARTIAL_BC: criteria (b) and (c) pass, (a) fails. Among those, dLPD median -19.03
  with SE 23.65, negative in 83%. Cause: fast-growing traces reach the oxygen floor
  partway through the window. With five CONTIGUOUS blocks, holding out the final
  block makes the exponential model EXTRAPOLATE beyond the crash, where it predicts
  physically impossible negative oxygen and is scored catastrophically, while the
  linear model degrades gracefully. This penalises exactly the traces with the
  clearest growth. It is an artefact of the CV geometry, not evidence about growth.

  OBSERVATION D (not repaired; NOT a defect of the classifier). Most remaining
  ambiguity sits at rt_true in (0, 0.5), where reason is PARTIAL_AB: (a) and (b)
  detect real curvature — correctly, since r > 0 by construction — while (c), the
  practical effect-size rule, fails. Reporting AMBIGUOUS there is the protocol
  behaving as designed on genuinely boundary inputs.

## V3 CHANGE (one, confined to criterion (a))

  CHANGE 3. Blocked predictive CV is made interpolative and physically bounded:
    (i)  Only the three INTERIOR blocks (2nd, 3rd, 4th of five) are held out. The
         training set therefore always brackets the held-out block in time, so no
         model is ever scored on extrapolation beyond the fitted range.
    (ii) Predictions are clipped to the physically attainable oxygen range
         [0, max(observed O2)] before the log predictive density is computed.
         Dissolved oxygen cannot be negative; scoring a model on impossible
         predicted values measures the extrapolation, not the fit.

  Both are justified by the measurement's support and the CV geometry, not by their
  effect on any acceptance number.

  UNCHANGED: every threshold, all three models, the respiration gate, criteria
  (b) and (c), R_FLOOR, RT_MIN, the four states, the G mapping.

## ACCEPTANCE CRITERIA: NOT MODIFIED

A1-A6 are carried over verbatim from v1.1. In particular A5 (pooled ambiguity
<= 0.25) is NOT relaxed, although OBSERVATION D indicates it is partly a property
of the chosen r-grid rather than of the classifier: a grid with substantial mass in
rt (0, 0.5) will produce correct-but-ambiguous calls there. If A5 still fails after
CHANGE 3, that fact is reported as-is, together with the ambiguity rate computed
separately among unambiguous truths (rt = 0 and rt >= 1), and the decision on how to
proceed is referred to the user rather than taken by relaxing the threshold.
