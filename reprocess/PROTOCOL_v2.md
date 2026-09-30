# PROTOCOL v2 — written after v1.1 FAILED simulation calibration
# Labels remain sealed. blind_key.csv has NOT been read.
# Parents: PROTOCOL_v1.md      f9aeed8f47317f5aae87bc9ca2ecb3dad61bf90f460a1e60f7be4c80ed8d3c2a
#          PROTOCOL_v1.1_AMEND 990320e3d468d8da52cec4bd91227c0965b6680197221111ff37da0caec8bf59

## V1.1 FAILURE RECORD (simulation, 1800 traces, sim_results.csv)

  A1 false growth at r=0    0.000 <= 0.05   PASS
  A2 power at rt >= 1       0.000 >= 0.80   FAIL
  A3 STATE2 at rt ~ 0.5     0.000 <= 0.50   PASS
  A4 STATE2 at rt >= 1      0.000 <= 0.10   PASS
  A5 ambiguous overall      0.967 <= 0.25   FAIL
  A6 STATE4 when dd > 2     0.000 <= 0.05   PASS

Diagnosis: reason code LOW_EFF_N on 1,726 of 1,800 traces. The classifier rejected
almost every trace before reaching any growth criterion, so power collapsed to zero
and nearly everything fell to AMBIGUOUS. Two distinct errors in v1:

  ERROR A (measurement). rho was the lag-1 autocorrelation of the LINEAR (M1)
  residuals. For any curved trace those residuals are dominated by model MISFIT,
  not noise, forcing rho -> 0.99 whatever the true noise structure. The observed
  median rho of 0.990 is therefore largely an artefact of fitting a straight line
  to a curve.

  ERROR B (conceptual, the more serious). eff_n = n(1-rho)/(1+rho) is the effective
  sample size for estimating a MEAN under AR(1). It is not the information content
  of a smooth, densely sampled trace about the parameters of a non-linear curve. A
  1-minute-resolution oxygen trace with 0.009 mg/L measurement noise carries a great
  deal of information about r and K even though adjacent points are near-perfectly
  correlated. Using this quantity as an admission gate was wrong in kind, not in
  calibration, and no adjustment of MIN_EFF_N repairs it.

## V2 CHANGES (two, both confined to the autocorrelation treatment)

  CHANGE 1. The MIN_EFF_N admission gate is REMOVED. eff_n is retained and reported
  for every well as a descriptor only. Justification: v1 already specifies two
  principled treatments of temporal autocorrelation — blocked predictive CV over
  five contiguous blocks (criterion a) and a parametric bootstrap with AR(1)-matched
  residuals (criterion b). Both handle correlation correctly and neither relies on
  eff_n. The gate was redundant as well as wrong.

  CHANGE 2. rho is estimated from the residuals of the BEST-FITTING model (M2 where
  it converges, M1 otherwise) rather than always from M1, so that model misfit is
  not mistaken for noise. rho is used only to generate AR(1)-matched bootstrap
  residuals and to report eff_n.

  UNCHANGED: usable-interval rules; decimation; the three models; the respiration
  gate (1.0 mg/L); all three growth criteria (a), (b), (c); R_FLOOR = 0.005 h^-1;
  RT_MIN = 0.5; the four states; the G = 0 / 1-2 / 3 mapping; all six acceptance
  criteria A1-A6, which are carried over UNMODIFIED and were NOT relaxed to secure
  a pass.

## RE-RUN REQUIREMENT

The full 1,800-trace simulation is rerun under v2 and must satisfy A1-A6 as
originally written. If it fails again, v2 is marked failed and a v3 is written
while labels remain sealed. Under no circumstances are acceptance thresholds
weakened to obtain a pass.

## V2.1 COMPUTATIONAL NOTE (appended before any v2 result was inspected)

Bootstrap refits start from the original point estimate (single start) rather than
the three-restart sweep used for the primary fit. A bootstrap replicate perturbs the
data only slightly, so its optimum lies near the original; the multi-start sweep is
needed for the primary fit and is retained there. This changes no threshold, no
model and no decision rule; it reduces bootstrap cost roughly threefold. r remains
bounded below by exactly 0 so the boundary stays attainable.
