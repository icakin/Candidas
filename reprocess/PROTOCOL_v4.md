# PROTOCOL v4 — written after v3 passed 5/6 but produced invalid real-data output
# Labels remain sealed. blind_key.csv has STILL not been read.

## V3 STATUS
Simulation: A1 PASS(0.000) A2 PASS(0.981) A3 PASS(0.000) A4 PASS(0.000)
            A5 FAIL(0.328 vs 0.25) A6 PASS(0.000)
Ambiguity among unambiguous truths (r=0 or rt>=1): 0.043. Ambiguity at rt in
(0,0.5), where the truth is genuinely borderline: 0.844. A5 therefore measures the
composition of the r-grid, not the classifier. NOT relaxed; see v3 sec. on deferral.

## TWO DEFECTS FOUND ON THE REAL TRACES (before any label was joined)

DEFECT E — the absolute oxygen floor almost never fires.
  O2_FLOOR was frozen at 0.5 mg/L on the assumption that this is the sensor's
  anoxic bound. The data contradict it: the minimum oxygen reached within the
  usable window has median 0.86 mg/L, and 77.5% of traces never fall to 0.5.
  Consequently the window ran to the 15 h cap for most wells, and a MEDIAN OF 49.9%
  of each window was post-depletion plateau. Fitting the exponential model across a
  window half of which is flat drives r -> 0: in the v3 real-data output, states 2
  and 3 have median rt = 0.000 and median dur_h = 14.99, against 6.29 h for state 1.
  The classifier was measuring the plateau, not the growth phase.

DEFECT F — the simulator did not reproduce that plateau.
  Simulated traces decayed exponentially to the target drawdown across the whole
  window and were clipped only at 0.02 mg/L, so they had no flat tail. The v3
  acceptance results therefore certify the classifier on traces that do not
  resemble the real ones. This invalidates the v3 calibration regardless of A5.

## V4 CHANGES

  CHANGE 4 (window end, repairs E). The usable interval ends at the first time
    O2 <= ymin + DEPLETION_FRAC * (O2_start - ymin),  DEPLETION_FRAC = 0.02
  i.e. when depletion is 98% complete, whatever the absolute oxygen value. The
  absolute 0.5 mg/L floor is retained only as a backstop, whichever comes first.
  Justification: DEPLETION_FRAC = 0.02 is the legacy pipeline's own trimming rule
  (2% above the series minimum) and therefore PREDATES this audit, the same class
  of justification as RT_MIN = 0.5.

  CHANGE 5 (simulator realism, repairs F). Simulated traces are clipped at a
  plateau floor drawn from the observed distribution of per-trace oxygen minima,
  so that simulated and real traces share the same post-depletion geometry. This
  changes no classifier threshold; it makes the calibration valid.

  UNCHANGED: equilibration start, decimation, the three models, the respiration
  gate, criteria (a) interior-block clipped CV, (b) R_FLOOR, (c) RT_MIN, the four
  states, the G mapping, and acceptance criteria A1-A6 verbatim.

## RERUN
Full 1,800-trace simulation under v4, then the 1,380 real wells. A5 is expected to
remain grid-dependent; the decision on A5 stays with the user and is not taken here.
