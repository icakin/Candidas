# FROZEN CLASSIFICATION PROTOCOL v1
# Automatic, blind re-derivation of oxygen-trace states
# Written BEFORE any new biological comparison. Freeze this file, then run.

## 0. DECLARATION OF PARTIAL UNBLINDING (limitation, stated honestly)

The analyst writing this protocol has already examined taxon-level patterns in this
dataset during the preceding audit, including which taxa lose growth at which
temperatures. Naive blinding is therefore NOT achievable and must not be claimed.

Mitigations, all of which are verifiable from this document:
  (a) Every threshold below is fixed either to a value that PREDATES this audit
      (the pipeline's own rt > 0.5 rule; the manuscript's 1.5 h equilibration
      exclusion) or to a physical/instrumental justification stated inline.
  (b) No threshold is chosen by reference to any taxon contrast.
  (c) The classifier code receives ONLY a random well ID and the trace. Taxon,
      isolate, temperature and replicate are withheld from the classifying
      function. Temperature is used nowhere in classification.
  (d) This file is hashed and timestamped before the classifier is run, and
      before biological labels are rejoined.
  (e) Sensitivity to neighbouring thresholds is reported, but the frozen values
      define the primary result regardless of the outcome.

## 1. DATA SOURCE

Input: results/tables/Oxygen_All_Long.csv  (2,033,640 rows; 1,380 series;
untrimmed, unsmoothed; columns Time [min], T, OTU, Replicate, Oxygen [mg/L]).
This is the earliest tabular form of the SensorDish exports in the repository.
Raw per-plate CSVs in data/ are the upstream source and are used only to verify
that Oxygen_All_Long preserves every timepoint (spot-check >= 10 series).

Nothing in the legacy pipeline is modified, deleted or overwritten. All outputs
of this protocol are written under reprocess/ only.

Blinding: a lookup table maps (T, OTU, Replicate) -> random WELLID (UUID4), written
to reprocess/blind_key.csv. The classifier reads reprocess/blind_traces.parquet
containing ONLY (WELLID, Time, Oxygen). The key is not read until section 7.

## 2. PRESPECIFIED USABLE INTERVAL  (no per-well hand selection)

  start_time   = max(90 min, argmax of Oxygen within the first 120 min)
                 Justification: 90 min is the manuscript's 1.5 h thermal-
                 equilibration exclusion. The reader's optical warm-up produces a
                 rising drift peaking at ~43-52 min (documented in the project's
                 own data dictionary); taking the later of the two removes both.
  oxygen_floor = 0.5 mg/L. The interval ends at the first time Oxygen <= 0.5,
                 because below this the sensor is at its anoxic bound and the
                 trace carries no rate information.
  max_duration = 15 h from start_time (caps late-phase / re-aeration artefacts).
  min_duration = 2 h.
  min_points   = 60 raw observations within the interval.
  min_eff_n    = 10 effectively independent observations (section 4).
  early anoxia = if the floor is reached before min_duration, the interval is
                 start_time -> anoxia time, and the trace is eligible only if it
                 still has >= min_points and >= min_eff_n; otherwise AMBIGUOUS
                 with reason code EARLY_ANOXIA_SHORT.

No smoothing is applied. Models are fitted to raw observations.

## 3. COMPETING TRACE MODELS  (identical interval, identical data)

  M0  constant      O2(t) = c
  M1  linear        O2(t) = c - b*t,                b >= 0
  M2  exponential   O2(t) = O2_0 + (K/r)(1 - exp(r*t)),   K > 0, r > 0
                    (the manuscript's mechanistic model)

M1 is the r -> 0 boundary of M2. The null r = 0 therefore lies ON the parameter
boundary, so the likelihood-ratio statistic is NOT chi-square distributed. The
comparison is made by blocked predictive CV (primary) and parametric bootstrap
(confirmatory); no chi-square p-value is computed.

## 4. AUTOCORRELATION AND SAMPLE SIZE

Residuals of M1 are strongly temporally autocorrelated; raw n (~1,400) grossly
overstates information. Procedure:

  rho      = lag-1 autocorrelation of M1 residuals.
  eff_n    = n * (1 - rho) / (1 + rho), reported for every well.
  PRIMARY  = blocked predictive cross-validation: the usable interval is split
             into 5 contiguous, equal-duration blocks. For each block, fit on the
             other four, predict the held-out block, score by mean log predictive
             density under a Gaussian with the training residual sd.
             dLPD(M2 - M1) is the primary statistic; its uncertainty is the SE
             across the 5 blocks.
  AICc     may be reported using eff_n in place of n, as a SECONDARY descriptor.
             It does not determine any state.

## 5. CLASSIFICATION RULES  (frozen)

Respiration gate:
  R : total oxygen drawdown across the usable interval >= 1.0 mg/L.
      Justification: instrument resolution and drift over a multi-hour window;
      calibrated against blank/sterile wells in section 6 if any exist.

Growth criteria, ALL evaluated only when R holds:
  (a) PREDICTION : blocked-CV dLPD(M2 - M1) > 0 with its across-block SE
                   excluding 0 (mean > 1 SE above 0).
  (b) IDENTIFIED : the 95% lower bound of r, from parametric bootstrap
                   (200 resamples with AR(1)-matched residuals), is > 0.
  (c) EFFECT SIZE: rt = r_hat * usable_duration_h >= 0.5.
                   Justification: this is the pipeline's OWN pre-existing
                   CURV_MIN_RT rule and predates this audit.

Let G = number of (a),(b),(c) satisfied.

  R false           -> STATE 4  NO DETECTABLE RESPIRATION
  R true,  G = 3    -> STATE 1  GROWTH + RESPIRATION
  R true,  G = 0    -> STATE 2  RESPIRATION, NO DETECTABLE GROWTH
  R true,  G = 1,2  -> STATE 3  AMBIGUOUS (conflicting evidence)
  interval invalid  -> STATE 3  AMBIGUOUS (reason code recorded)

Conflicting criteria NEVER produce a biological state. A positive fitted r alone
never produces STATE 1.

Reason codes recorded for every well: OK, SHORT_WINDOW, FEW_POINTS, LOW_EFF_N,
EARLY_ANOXIA_SHORT, NO_SIGNAL, FIT_FAILED, PARTIAL_EVIDENCE_A, _B, _C.

## 6. CALIBRATION BEFORE UNBLINDING

Simulation: traces generated across the observed ranges of O2_0, K, duration and
noise, with AR(1) residuals matched to observed rho, at r in
{0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8} h^-1. >= 200 traces per r.
All simulated traces pass through the identical classifier.
Report: false-positive STATE 1 rate at r = 0; power to call STATE 1 by true r;
the practical detection boundary (smallest r called STATE 1 in >= 80% of traces);
ambiguous rate by r.
Blanks/sterile controls, if present in the raw exports, calibrate the R gate.

## 7. FREEZE AND UNBLIND

Frozen before unblinding: this file's hash, the classifier code's hash, all
thresholds, the simulation performance table, and reprocess/blind_states.csv
(WELLID -> state + reason code + all statistics).

Only then is blind_key.csv joined. Post-unblinding outputs, all prespecified:
  - confusion matrix, automatic state vs legacy manual include/exclude;
  - state counts by taxon, isolate and temperature, over all 1,200 non-glabrata;
  - ambiguous counts by cell;
  - whether relatives still enter STATE 2 at 40 C, and whether that conclusion
    changes if all AMBIGUOUS wells are reassigned to STATE 1 (worst case) ;
  - isolate-level (not only well-level) support, counting isolates not wells.

Figures 1-3 are NOT refitted and no biological contrast is interpreted until the
above is complete.
