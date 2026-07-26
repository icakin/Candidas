# runs/ — MANIFEST

`runs/` holds **non-canonical, per-prompt output trees**. Each is the complete
record of one analysis run, kept because a report cites it.

`results/` is the **canonical, shipped** tree. Nothing here replaces it, and
nothing here is read by the pipeline unless a `CANDIDAS_*` environment variable
points at it explicitly.

> **Every tree is auditable from this file alone.** Each entry records what
> produced it, when, on what environment, what it is evidence for, and a
> checksum index of its tables — so a tree can be pruned later and its integrity
> still verified against what was reported.

**Environment for every run below** — full detail in
[`env/versions.json`](../env/versions.json):

| | |
|---|---|
| machine | macOS 26.2, Apple silicon (arm64), 16 cores |
| R | 4.5.2, `renv.lock` (170 packages, CRAN + Bioconductor 3.22) |
| Stan | rstan 2.32.7 / StanHeaders 2.32.10 / Stan 2.32.2 (brms default backend) |
| Python | 3.9.6, `cauris_etcgem/requirements.lock` |
| etc-GEM engine | `etcgem @ git+…/etcGEMs@7d32383d902184c085662031b802406353e7b083` |
| LP solver | GLPK 5.0 via swiglpk 5.0.13 |
| Quarto | 1.8.27, TinyTeX v2026.02 |

---

## `C1_reproduction/`

| | |
|---|---|
| **produced by** | `prompts/C1_environment_and_end_to_end_reproduction_prompt.md` |
| **command** | `bash scripts/run_c1.sh` — i.e. `bash scripts/run_all.sh` with `CANDIDAS_RESULTS=<tree>`, `ETCGEM_OUT_SUFFIX=_C1`, `CANDIDAS_SUPP_DATA=<etc-GEM C1 tree>`, `CANDIDAS_EXPRESSION_OUT=<tree>/expression` |
| **when** | 2026-07-26, 09:47:44 → 10:36:34 (48 m 50 s: etc-GEM 35 m 29 s, R 12 m 37 s, render 44 s) |
| **commit** | `ab1686d` — *C1: rerun the R pipeline into results_C1 + render the manuscript* |
| **N0 treatment** | shipped (`N0_BACKPROJECT = TRUE`), no seed |
| **logs** | `logs/run_all_20260726_094744.log`, `logs/{etcgem,r_pipeline,quarto}_20260726_094744.log` |

**Evidence for** — [`reports/REPRODUCTION.md`](../reports/REPRODUCTION.md), all
sections. Specifically:

* §1a, the file-by-file comparison of `results/` against this tree — the
  finding that everything through stage 07 reproduces IDENTICAL or NUMERICALLY
  EQUIVALENT (worst 2.36e-10), and everything downstream of stage 08 is
  STOCHASTIC-BUT-CONSISTENT;
* §2, the headline table of load-bearing numbers;
* §8, run-to-run determinism (39 of 40 tables byte-identical across two
  complete runs).

It is **also the ARM 1 baseline for C2**: the reproducibility floor in
`reports/N0_SENSITIVITY.md` §3.0 is this tree against `C2_arm_current/`, and
every PART A/B diagnostic in that report reads its
`tables/derived_N0_R_results_with_carbon.csv`,
`Oxygen_Trimmed_Series_Metadata.csv`, `Oxygen_All_Long.csv` and
`sharpe_schoolfield_growth_fgC_h_coefs.csv`.

| contents | tracked? |
|---|---|
| `tables/` — 37 CSVs (+3 large oxygen intermediates, gitignored by name) | **yes** |
| `expression/` — the GEO differential-expression output of stage 11 | **yes** |
| `figures/manuscript/` — the 7 publication figures | **yes** |
| `figures/` — 71 further per-series diagnostic figures | no (regenerable) |
| `rds/` — 5 cached brms fits, 29 MB | no (regenerable) |
| `manuscript_build/` — staging render against these figures, 20 MB | no (regenerable) |
| `CHECKSUMS_tables.txt` — SHA-256 of all 37 tables | **yes** |

---

## `C2_arm_current/`, `C2_arm_ramp/`, `C2_arm_nobp/`

| | |
|---|---|
| **produced by** | `prompts/C2_n0_backprojection_decision_prompt.md` |
| **command** | `bash scripts/run_c2_arms.sh` |
| **when** | 2026-07-26, 11:18 → 11:57 (current 10 m 33 s, ramp 10 m 23 s, nobp 11 m 07 s) |
| **commit** | `ada23b8` — *C2: run three N0 arms end to end* |
| **logs** | `logs/c2_arm_{current,ramp,nobp}_20260726_111821.log` |

Each is a complete re-run of stages 07 → 08 → 09 → 10 → 12 → 13 with a full
Bayesian refit. **The only thing that differs between them is the $N_0$
treatment**:

| tree | `CANDIDAS_N0_MODE` | $N_0$ | status |
|---|---|---|---|
| `C2_arm_current/` | `current` | $N_{\text{inoc}}e^{r\delta}$ | the shipped treatment |
| `C2_arm_ramp/` | `ramp` | $N_{\text{inoc}}e^{r\delta}f_{\text{ramp}}$ | a **bound**, not an estimator |
| `C2_arm_nobp/` | `nobp` | $N_{\text{inoc}}$ | a **lower bound on $N_0$**, not a candidate |

Held identical across all three: `CANDIDAS_SEED=20260726`; the committed trim
windows, exclusions, inoculum and cell-carbon constants; every model formula,
prior and the sampler seed 1234; the etc-GEM capacity values (all three read
`cauris_etcgem/runs/C1/supp_data/`); and the stage 02/03/06 outputs, which are
$N_0$-independent and were carried over byte-for-byte from `C1_reproduction/`.
The ARM 2 multipliers come from
[`reports/tools/c2_ramp_n0_factors.csv`](../reports/tools/c2_ramp_n0_factors.csv).

**Evidence for** — [`reports/N0_SENSITIVITY.md`](../reports/N0_SENSITIVITY.md)
and [`reports/n0/n0_report.pdf`](../reports/n0/n0_report.pdf).
Specifically:

* the invariance check — every growth-side quantity differs by **exactly
  0.000e+00** between arms, which is what shows the switch does not leak;
* the $E_R$, $\Delta E$, $T_{\text{opt}}$(CUE) and carbon-tax tables;
* the ten pairwise fever contrasts, and the finding that Clade III vs Clade IV
  resolves in `current` and `ramp` but not in `nobp`;
* the capacity correlations, and the collapse of $r$(capacity, respiration)
  from $-0.620$ to $-0.196$ ($p = 0.54$) in `nobp`.

| contents (each tree) | tracked? |
|---|---|
| `tables/` — 37 CSVs | **yes** |
| `figures/manuscript/` — the publication figures for that arm | **yes** |
| `figures/` — 69 further per-series diagnostic figures | no (regenerable) |
| `rds/` — 5 cached brms fits, 29 MB | no (regenerable) |
| `CHECKSUMS_tables.txt` — SHA-256 of all 37 tables | **yes** |

---

## `cauris_etcgem/runs/C1/`

Inside the submodule, because that is where the etc-GEM pipeline writes.

| | |
|---|---|
| **produced by** | `prompts/C1_environment_and_end_to_end_reproduction_prompt.md`, PART C |
| **command** | `generate_model_data.py all --out-suffix _C1` (via `scripts/run_c1.sh`) |
| **when** | 2026-07-26, 09:47:47 → 10:23:16 (35 m 29 s) |
| **commit** | submodule `ead1e20` — *C1: rerun the etc-GEM pipeline into supp_data_C1 (incl. stage_boot)* |
| **stage timings** | build 2.5 s · calibrate 16 m 00 s · curves/fit < 1 s · ident 61.7 s · boot 92.6 s · apriori 8.6 s · bayes 16 m 26 s |
| **log** | `logs/etcgem_20260726_094744.log` |

| contents | |
|---|---|
| `supp_data/` | 13 CSVs + `bayes_chain_raw.npz` (10 MB: the full un-thinned emcee chains, walker- and step-indexed, with per-walker acceptance fractions and autocorrelation times) |
| `figure3_data/` | 4 CSVs |
| `_percladefit.npy` | the differential-evolution fit cache — **gitignored**, regenerable, and deliberately so: committing it would let a stale calibration be silently inherited, which is one of the failure modes C1 documented |
| `CHECKSUMS.txt` | SHA-256 of all 19 files |

**Evidence for** — `reports/REPRODUCTION.md` §1b, §1c and §2d: that the etc-GEM
outputs do **not** reproduce, and specifically that `kcat_scale`
0.614 / 0.410 / 0.622 / 0.710 regenerates as 0.833 / 0.561 / 0.851 / 0.945, a
rescale of 1.33–1.37× against an anchor ratio of 1.0890/0.7955 = 1.3690.

It is **also the fixed capacity axis for all three C2 arms**, which is why it
must not be regenerated: C1 established the etc-GEM is not reproducible (emcee
is unseeded; differential evolution is seeded but its objective is an LP solve),
so letting capacity vary would confound the arm comparison.

---

## What must never be deleted

| path | why |
|---|---|
| `results/` | the canonical shipped tree; the reproduction comparison is against it |
| `data/` | the raw instrument exports |
| `cauris_etcgem/strains/eci_cauris/outputs/supp_data/` | the committed etc-GEM outputs |
| `cauris_etcgem/strains/eci_cauris/outputs/supp_data_original_backup/` | **the mixture of these two directories is itself the evidence** for the provenance finding in `reports/REPRODUCTION.md` §4.3 — four of thirteen files are byte-identical between them and the rest are not |
| `cauris_etcgem/strains/eci_cauris/outputs/figure3_data/` | evidence for §4.4 |
| `manuscript/` | the write-up |

---

## Regenerating any tree

```bash
# C1_reproduction/ and cauris_etcgem/runs/C1/          (~49 min)
bash scripts/run_c1.sh

# C2_arm_current/, C2_arm_ramp/, C2_arm_nobp/          (~32 min)
bash scripts/run_c2_arms.sh                 # or: bash scripts/run_c2_arms.sh ramp
```

Both scripts are non-destructive: they write only to their own tree and never
to `results/`, `data/` or the committed etc-GEM outputs. Verify with:

```bash
shasum -a 256 -c env/baseline_checksums_results.txt
shasum -a 256 -c env/baseline_checksums_etcgem.txt
shasum -a 256 -c env/baseline_checksums_data.txt
```

Verify a run tree against what was reported:

```bash
cd runs/<tree>/tables && shasum -a 256 -c ../CHECKSUMS_tables.txt
```
