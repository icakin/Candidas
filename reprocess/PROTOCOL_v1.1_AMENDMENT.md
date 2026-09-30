# PROTOCOL v1.1 — AMENDMENT TO v1 (made while still blind; no labels inspected)
# Parent: PROTOCOL_v1.md  sha256 f9aeed8f47317f5aae87bc9ca2ecb3dad61bf90f460a1e60f7be4c80ed8d3c2a
# Reason: a pre-run completeness check found two defects in v1. Both are repaired
# here BEFORE any classifier is executed and before any label is joined.

## DEFECT 1 — criterion (b) was vacuous.  [REPAIRED]

v1 required "the 95% lower bound of r, from parametric bootstrap, is > 0".
In v1 the optimiser constrained r to a strictly positive lower bound, so every
bootstrap replicate necessarily returned r > 0 and the criterion could never fail.
It carried no information.

Repair, two parts:
  (i)  The optimiser lower bound for r is set to EXACTLY 0, so r = 0 is attainable
       and bootstrap replicates can genuinely pile up at the boundary.
  (ii) Criterion (b) becomes: the 5th percentile of the bootstrap distribution of
       r exceeds R_FLOOR = 0.005 h^-1.
       Justification for R_FLOOR: one twentieth of the slowest median growth rate
       reported anywhere in the legacy pipeline (~0.1 h^-1), and below any rate
       that could produce one e-fold within the longest usable window. It is a
       numerical-noise floor, not a biological threshold, and it is deliberately
       an order of magnitude below criterion (c).

## DEFECT 2 — no simulation acceptance criteria were specified.  [REPAIRED]

v1 stated what the simulation must report but never what performance would make
the protocol fit to use. Acceptance criteria are now fixed, in advance:

  A1  False growth: STATE 1 called for <= 5% of simulated traces with true r = 0.
  A2  Power: STATE 1 called for >= 80% of traces with true rt >= 1.0
      (i.e. one e-fold across the usable window, twice criterion (c)).
  A3  Boundary behaviour: at true rt = 0.5, STATE 2 ("respiration, no detectable
      growth") is called for <= 50% of traces. At the criterion boundary the
      classifier must not confidently assert absence of growth.
  A4  Misclassification: STATE 2 called for <= 10% of traces with true rt >= 1.0.
  A5  Ambiguity: STATE 3 assigned to <= 25% of all simulated traces pooled across
      the r grid.
  A6  Respiration gate: STATE 4 ("no detectable respiration") called for <= 5% of
      traces whose simulated drawdown exceeds 2.0 mg/L.

  ALL SIX must hold. If any fails, protocol v1.1 is marked FAILED, the failure is
  documented, a justified v2 is written and hashed while labels remain hidden, and
  the simulation is rerun. Under no circumstances are labels joined to a failed
  classifier's output.

## CLARIFICATIONS (v1 text was adequate but is restated numerically)

  (1) "Better blocked prediction" = mean dLPD(M2 - M1) across the 5 held-out
      blocks > 1.0 x its standard error across those blocks, and > 0 in absolute
      terms. Both conditions.
  (3) Respiration-only vs inert is decided solely by the R gate: drawdown over the
      usable interval >= 1.0 mg/L is respiring; < 1.0 mg/L is STATE 4.
  (4) Disagreement among (a), (b), (c): G = count satisfied; G = 3 -> STATE 1,
      G = 0 -> STATE 2, G in {1,2} -> STATE 3 AMBIGUOUS with the passing criteria
      recorded as reason codes. No tie-breaking, no ordering of criteria.

## COMPUTATIONAL NOTE (does not alter any threshold)

Traces carry ~1,400 points with strong AR(1) structure. For tractability every
trace is uniformly decimated to at most 300 points within the usable interval
BEFORE fitting. Decimation is uniform in time, applied identically to every well
and to every simulated trace, and is therefore not a selection rule. All reported
eff_n are computed post-decimation.
