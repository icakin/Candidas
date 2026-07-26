# HANDOVER — Candidas, after C1, C2 and C2b

For whoever picks this up next, including its author returning cold. Everything
links out rather than restating.

---

## What is solid

**The pipeline runs.** A clone plus `Rscript scripts/00_install.R` and
`bash scripts/run_all.sh` regenerates everything from the committed raw data in
**48 min 50 s**, unattended, no browser, no prompts. Before C1 it could not run
at all: `run_all.R` sourced `06_inoculation.R`, which called `runApp()` under
`Rscript`, so the master runner blocked forever on a Shiny server nobody could
see.

**The respirometry half reproduces.** Every table through stage 07 comes back
IDENTICAL or NUMERICALLY EQUIVALENT — worst disagreement **2.4 × 10⁻¹⁰**. Every
load-bearing Bayesian number returns within **0.4 % relative**, inside its own
credible interval. The respiration Arrhenius fit reproduces *exactly*: every
posterior mean, sd, quantile and R̂ differs by 0.0.

**The N₀ treatment does not touch growth.** Across all three N₀ arms, `r`, `K`,
absolute and mass-specific growth, the growth posterior, the per-isolate growth
coefficients, E_G and the growth optimum all differ by **exactly 0.000 × 10⁰**.
The switch does not leak.

**The headline result holds, and strengthens.** *Every taxon's carbon-economy
optimum lies below body temperature* with posterior probability **1.0000 in all
fifteen taxon-arm combinations**. It is robust to a treatment that moves E_R by
up to 0.32 eV. The paper currently understates this.

---

## Status

**C1 — reproduction.** [`reports/REPRODUCTION.md`](reports/REPRODUCTION.md).
Pinned R/Python/Stan/Quarto (`renv.lock`, `requirements.lock` with the etc-GEM
engine at a commit), made the pipeline headless, re-ran everything into a
parallel tree. Found and confirmed all seven predicted internal inconsistencies,
and identified the cause of the etc-GEM failure: the committed
`calibration_r2.csv` was built with the **previous anchor**
(`ANCHOR_MU = 1.1 @ 313.15 K`), and 1.089/0.7955 = **1.3690** is exactly the
observed capacity rescale.

**C2 — the N₀ back-calculation.** [`reports/N0_SENSITIVITY.md`](reports/N0_SENSITIVITY.md)
and the rendered report [`reports/n0/n0_report.pdf`](reports/n0/n0_report.pdf)
(17 pp). `respiration = K/N₀` exactly, so
`log R = log K − log N_inoc − r·δ` (verified to sd 1.3 × 10⁻¹⁴). δ is a
*detection* delay, not a lag. Three arms re-run end to end.

**C2b — report and tidy.** The run trees are consolidated under
[`runs/`](runs/MANIFEST.md) with a manifest giving each one's producer, command,
environment, purpose and table checksums. 654 MB of provably regenerable,
uncited bulk pruned; every table kept.

---

## Manuscript status, claim by claim

| claim | verdict | the number |
|---|---|---|
| Growth turns over; respiration does not | **SURVIVES** | LOO prefers Arrhenius over Sharpe-Schoolfield for respiration in every arm; unchanged by N₀ |
| Every economic optimum lies below 37 °C | **SURVIVES** — strengthen it | P = **1.0000**, all 5 taxa × 3 arms |
| All four *C. auris* clades pay credibly less at fever than *C. parapsilosis* | **SURVIVES** | 4/4 contrasts resolve in **all three** arms (12/12) |
| Growth is ~2× more temperature-sensitive than respiration | **SURVIVES as direction, WEAKENS as "2×"** | E_G/E_R median 2.36 (current), 2.94 (ramp), **1.58** (lower bound). ΔE still excludes zero in all five taxa everywhere |
| **Clade IV is credibly the cheapest at fever** | **FAILS** | rests on one contrast. Clade III vs IV: 1.198 [1.031, 1.394] current, 1.176 [1.029, 1.352] ramp, **1.100 [0.920, 1.319]** lower bound — CrI spans 1 |
| Spearman ρ = −0.90 [−1.00, −0.30] | **WEAKENS — two separate problems** | (i) it does not match its own results: the committed fits give **−1.00 [−1.00, −0.70]**, twice, and so do all three arms. (ii) n = 5, and both axes (growth T_opt and fever cost) are functions of the *same* four posterior shape parameters, so it is not an independent correlation. Not an N₀ artefact — it is invariant across arms |
| Fig 4e, capacity vs per-cell respiration r = −0.64 | **FAILS under the lower bound** | −0.620 → −0.603 → **−0.196, p = 0.54**. ARM 3 respiration contains no `r`; capacity is fitted to growth. The correlation exists only through `e^(r·δ)`. Fig 4d (−0.85) is stable and does not have this problem |

> The prompt for this handover quoted the failing contrast as "1.00 [0.920,
> 1.319]". The measured value is **1.100**; the interval is right. Using the
> measured one.

---

## What is open

**1. The N₀ back-calculation — a method-level question, not a Candidas one.**
The same construction is in `icakin/OxygenModel` (`CUE.R:222`). φ, the true
growth over δ, is **not identified** by the oxygen record: it enters only as a
multiplicative constant on N₀, which is not a parameter of the fitted curve, so
the profile likelihood in φ is exactly flat. The three arms are three *priors*,
not three estimates ([`n0_report.pdf`](reports/n0/n0_report.pdf), Appendix).
*Next action:* it belongs in OxygenModel, where endpoint counts exist. **D2
there already found forward and backward N₀ differ by a median factor of 57**,
growing with `r` — so this is not a small correction. See
`../OxygenModel/HANDOVER.md`.

**2. `T_internal` — recorded and stripped.** All 72 raw PreSens exports carry a
measured temperature trace; `01_convert_xlsx.R` discards it, along with `Tm`.
The plate sits **4.75 K** below its settled temperature ten minutes into a 44 °C
run and is still 0.40 K low at 84 min, while the optode compensated at a
**constant** `Tm` throughout — so reported O₂ *concentrations* are compensated at
a temperature wrong by up to ~4.8 K over exactly the pre-window interval.
*Next action:* recompute concentrations from phase and the measured temperature.
Until then the fit-window contamination bound rests on a transient that is
itself partly a compensation artefact. **This also unblocks OxygenModel's open
item 4**, which is waiting on exactly these `.xlsx` files —
they are in [`data/`](data/), and
[`reports/tools/c2_extract_Tinternal.R`](reports/tools/c2_extract_Tinternal.R)
already extracts the trace.

**3. The 88 manual exclusions and the hot end.** 43 % are flagged by no
automatic criterion, they are directional, and *C. parapsilosis* has no data
above 40 °C yet is tabulated and plotted to 44 °C.
*Next action:* C3 — document each exclusion, model non-respiring wells as zeros
rather than as missing, stop extrapolating.

**4. The etc-GEM half — nothing absolute reproduces.** Capacity regenerates
1.33–1.37× higher; `outputs/supp_data/` is a **mixture of two builds** (four of
thirteen files byte-identical to the backup, the rest not); and
`figure3_data/model_curves.csv` uses a temperature grid the current code cannot
emit (32 points at 1 °C over 18–49 °C vs 53 at 0.5 °C over 18–44 °C).
*Next action:* deferred — see Decisions.

---

## The submodule problem

`cauris_etcgem` is **private**, and [`.gitmodules`](.gitmodules) points at
`https://github.com/icakin/cauris_etcgem.git` while the working remote is
`git@github.com:icakin/cauris_etcgem.git`. An unauthenticated
`git clone --recursive` therefore **fails**, and everything the etc-GEM figures
depend on is unreachable to anyone outside the group. This is the reproducibility
blocker C1 documented.

Two things would fix it — **neither has been done**, because changing the URL
without knowing which way it will go would make this note wrong:

* make `cauris_etcgem` public and the https URL works as declared; or
* vendor the etc-GEM outputs into this repository and drop the submodule.

---

## Decisions already taken

Not up for relitigation.

* **The back-calculation is a design choice forced by the sealed vial.** Cell
  density cannot be measured at the window start without destroying the
  closed-system measurement. The question is only *which* back-calculation.
* **ARM 3 (`N₀ = N_inoc`) is a lower bound, not a candidate.** δ is a detection
  delay; the cells were growing through it. It is retained as a *diagnostic*,
  because it is the only arm that removes the `r` dependence, so a result that
  moves between ARM 1 and ARM 3 is carried by the coupling rather than by the
  measurement.
* **`N0_BACKPROJECT` stays `TRUE`.** C2 changed no default; `CANDIDAS_N0_MODE`
  selects an arm and reproduces the shipped pipeline exactly when unset.
* **The etc-GEM decision is deferred** until `T_internal` is resolved and C3 has
  run. Whether it stays in the paper is not a C1/C2 call.
* **The estimator work moves to the OxygenModel repo**, where the endpoint
  biomass exists. There is **none** in this repository — every column of all 72
  raw exports was checked.

---

## How to resume

```bash
# install (10-40 min cold; verifies Stan by compiling a toy model)
Rscript scripts/00_install.R
cd cauris_etcgem && python3 -m venv .venv \
  && .venv/bin/python -m pip install -r requirements.lock && cd ..

# run everything into the canonical tree            (~49 min)
bash scripts/run_all.sh

# or non-destructively, into runs/
bash scripts/run_c1.sh          # full reproduction
bash scripts/run_c2_arms.sh     # the three N0 arms  (~32 min)
```

**Read in this order:** [`reports/n0/n0_report.pdf`](reports/n0/n0_report.pdf)
(the argument, 17 pp) → [`reports/REPRODUCTION.md`](reports/REPRODUCTION.md)
(what does and does not reproduce) →
[`reports/RUNBOOK.md`](reports/RUNBOOK.md) (what to check at each stage) →
[`runs/MANIFEST.md`](runs/MANIFEST.md) (which tree is evidence for what).

**Watch for:** the model build printing `mu* = 1.089` instead of `0.7955` — that
means the old anchor, and your numbers will match the committed
`calibration_r2.csv` rather than the current code. It is the single biggest
provenance trap here (`REPRODUCTION.md` §4.2).

[`prompts/README.md`](prompts/README.md) lists the run order. **C1, C2 and C2b
have been run. C3–C8 are described there but NOT YET WRITTEN** — the summaries
are a specification, not a prompt you can paste.

**Sister project:** `../OxygenModel` (`icakin/OxygenModel`) and its
`HANDOVER.md`. The two now share one open question — the N₀ back-calculation —
and each holds part of the answer: OxygenModel has the flow-cytometry endpoint
counts, Candidas has the raw `.xlsx` exports with `T_internal`.
