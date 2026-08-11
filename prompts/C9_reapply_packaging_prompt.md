# Claude Code prompt — C9: reapply the C1 packaging onto Ilgaz's main, against the new script numbering (autonomous)

Run from the project root (`.../Candidas TPC/Candidas`). Packaging and execution only. It changes NO
analysis decision, NO constant, NO script logic and NO number. The purpose is a small, reviewable pull
request that gets the pinned environment and the headless runner onto the current main, where four
new scripts and a renumbering have landed since our branch was cut.

NOTE TO USER: launch in an auto-approving mode. Needs network and credentials for
`git@github.com:icakin/Candidas.git`. The pipeline takes ~49 min, so budget most of that for the
verification run.

CONTEXT. `gyd/c1-c2-review` carries C1, C2 and C2b (23 commits). `origin/main` has since advanced by
seven commits of Ilgaz's own work: a manuscript v3 rebuild, four new analysis scripts, a renumbering,
vendored etc-GEM outputs, and preservation of `Tm`/`T_internal` in `01_convert_xlsx.R`.

DO NOT CHERRY-PICK. The numbering has shifted from 08 onward, so C1's `--file=` fixes — which went
into "02, 03, 07-10, 12" under the old scheme — would silently patch the WRONG scripts. Reapply by
hand against his layout.

    OURS (14)                        HIS (18, current main)
    01_convert_xlsx.R          ->    01_convert_xlsx.R          (he has since edited this)
    02_longdata.R              ->    02_longdata.R
    03_trimming.R              ->    03_trimming.R
    04_trim_selector.R         ->    04_trim_selector.R         (Shiny app)
    05_cell_sizes.R            ->    05_cell_sizes.R            (Shiny app)
    06_inoculation.R           ->    06_inoculation.R           (Shiny app)
    07_oxygen_fits.R           ->    07_oxygen_fits.R
    08_bayesian_models.R       ->    09_bayesian_models.R
    09_bayesian_plots.R        ->    11_bayesian_plots.R
    10_carbon_tax.R            ->    12_carbon_tax.R
    11_capacity_expression.R   ->    13_capacity_expression.R
    12_main_figures.R          ->    14_main_figures.R
    13_supplementary_figures.R ->    16_supplementary_figures.R
    14_schematic.py            ->    17_schematic.py
    (none)                     ->    08_temperature_equilibration_sensitivity.R   NEW
    (none)                     ->    10_n0_term_test.R                            NEW
    (none)                     ->    15_uncertainty_bands.R                       NEW
    (none)                     ->    18_n0_treatment_panel.R                      NEW

Treat that table as a starting point, not gospel: verify it against `origin/main` before relying on
it, and report any discrepancy.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: on `gyd/c1-c2-review` —
`SETUP.md`, `env/versions.json`, `env/baseline_checksums_*.txt`, `.Rprofile`, `.gitignore`,
`scripts/00_install.R`, `scripts/run_all.sh`, `scripts/run_all.R`, and the headless guards and
`--file=` resolution in the numbered scripts. Then on `origin/main` — `scripts/config.R`,
`scripts/run_all.R`, all eighteen numbered scripts, and `manuscript/draft/`. Understand both sides
before writing anything.

PART A - branch and local hygiene
- `git fetch origin`. Local `main` currently points at our handover commit rather than his; reset it:
  `git branch -f main origin/main`. Report that you have done so.
- Branch off his current head: `git switch -c gyd/packaging origin/main`. Print the base commit.
- Verify the script mapping table above against `origin/main` and report any difference.

PART B - reapply, adapted to his layout
Bring across ONLY these, rewritten for his numbering:
- `renv.lock` and the renv activation (`.Rprofile`, `renv/activate.R`, `renv/settings.json`). Verify
  the lockfile covers every package his EIGHTEEN scripts load, including anything new since we
  branched (the four new scripts may pull in packages our lockfile does not have — check
  `08_temperature_equilibration_sensitivity.R`, `10_n0_term_test.R`, `15_uncertainty_bands.R` and
  `18_n0_treatment_panel.R` specifically). Add what is missing; remove nothing.
- `env/versions.json` — REGENERATE on this machine rather than copying ours.
- `env/baseline_checksums_*.txt` — regenerate against HIS current `data/`, `results/` and etc-GEM
  outputs, since those have changed.
- `scripts/00_install.R` — idempotent installer, updated for the packages above.
- `scripts/run_all.sh` — `set -euo pipefail`, per-stage timing, tee'd log to `logs/`. It must call HIS
  script names and numbers, and must cover every non-interactive stage including the four new ones.
  Decide and state where 08, 10, 15 and 18 belong in the order; if any is a diagnostic rather than a
  pipeline stage, say so and leave it out with a comment explaining why.
- Headless guards on the three Shiny scripts (`04_trim_selector.R`, `05_cell_sizes.R`,
  `06_inoculation.R`) so they are sourceable without launching a browser. Reapply to HIS versions.
  Their committed outputs remain inputs and must not be regenerated.
- The `--file=` script-directory resolution, including the `~+~` un-mangling for paths containing
  spaces. Apply to every script that needs it under HIS numbering — work out which those are rather
  than copying our list.
- `SETUP.md` — rewritten for eighteen scripts and the current runtime.
- `.gitignore` — bring across our rules, and ADD `manuscript/.~lock.*#` (four LibreOffice lock files
  are currently tracked on his main). Keep LaTeX artefact rules scoped to the report directories; do
  NOT add a global `*.log`.

EXPLICITLY EXCLUDED — none of these may appear in the diff:
  reports/ · prompts/ · runs/ · HANDOVER.md · run_c1.sh · run_c2_arms.sh ·
  any results/ or data/ file · any change to config.R constants ·
  any change to model, fitting, figure or manuscript content ·
  the cauris_etcgem submodule pointer

PART C - verify, and this is the part that matters
- `bash scripts/run_all.sh` end to end, unattended. Report total and per-stage wall time. Every stage
  must run or be documented as intentionally skipped, with a reason.
- CRITICAL: confirm the run reproduces HIS committed `results/` unchanged. Compare every regenerated
  table against the committed copy and report max absolute and max relative difference per file. If
  ANY table moves beyond floating-point noise, STOP and report rather than committing — that would
  mean a packaging change altered behaviour, which is exactly what this PR must not do.
- Confirm his four NEW scripts run under the pinned environment. They have never been run by anyone
  but him; if one fails, report it as a finding with the error rather than fixing it.
- Confirm the manuscript still builds: `cd manuscript/draft && quarto render`.

PART D - the manuscript figure linkage: diagnose fully, and specify the fix

`manuscript/draft/` has been rebuilt: the hand-written `_body_main.qmd`, `_body_supp.qmd`,
`manuscript.qmd` and `supplementary.qmd` are gone, replaced by `v3_quarto/` containing extracted
bitmaps at `media/media/image1..10.png`. Previously the manuscript sourced figures from
`results/figures/manuscript/`, so regenerating the pipeline updated the paper.

THE PROJECT POSITION, which the PR body must state plainly rather than raise as an option: every
figure in the manuscript MUST resolve to a file produced by the pipeline, referenced by path.
Embedded or extracted bitmaps are not an acceptable durable arrangement, whatever the short-term
convenience of a Word round-trip. Two things break when the link is severed:

  * a regenerated figure never reaches the paper, so the manuscript can silently disagree with the
    analysis it reports - and C1 exists precisely because that had already happened elsewhere in this
    project;
  * a reader or referee cannot trace a figure back to the code that made it, which is the property
    the whole C-series has been buying.

Do the diagnosis here, in full, so the remediation is a small specified job rather than an
investigation:

- Report whether anything in `manuscript/draft/` still references `results/figures/`, and list any
  figure that is still a live input.
- INVENTORY every embedded image: for each of `media/media/image1..10.png`, identify the pipeline
  figure it corresponds to (match by dimensions, content and file size against
  `results/figures/manuscript/` and `results/figures/`), and give the producing script for each.
- FLAG ANY EMBEDDED IMAGE WITH NO PIPELINE COUNTERPART. Those are figures with no producing code, and
  they are the serious case - the same failure this project has already hit once. Name them
  explicitly and say what they appear to show.
- Write the remediation spec to `reports/MANUSCRIPT_FIGURE_LINKAGE.md`: the image-to-source mapping,
  the figures needing a producing script, and what the `.qmd` needs in order to reference
  `results/figures/manuscript/` by path again. Concrete enough that restoring it is a follow-up
  prompt, not a research task.

Do NOT change any manuscript file in this PR. The packaging change must stay small and obviously
safe; the linkage fix is the next prompt and needs its own verification that the rendered output is
unchanged. But state the requirement, not a suggestion, in both the PR body and the spec.

PART E - the pull request
- Push `gyd/packaging` and open a PR against `origin/main` with `gh`. Title: "Packaging: pinned
  environment, headless guards, complete runner". Body: what is included; that it changes no analysis
  decision, constant or number; the verification that his committed results regenerate unchanged; the
  runtime; whether his four new scripts ran; and the PART D observation.
- Do NOT push to `main`. Do NOT force-push. Verify `origin/main` unchanged at the end.

VERIFY (report all)
1. Local `main` reset to `origin/main`; base commit branched from; the script mapping verified or
   corrected.
2. Every file added or modified with a one-line reason; confirmation the excluded list is absent from
   `git diff origin/main...gyd/packaging --stat`.
3. `bash scripts/run_all.sh` completes unattended; total and per-stage time; every stage accounted for.
4. His committed `results/` regenerate unchanged — per-file max absolute and relative difference.
5. The four new scripts: ran / failed, with errors verbatim if they failed.
6. `renv.lock` covers all eighteen scripts' dependencies; what was added.
7. Manuscript renders. The full PART D inventory: which images map to which pipeline figure and
   producing script; any embedded image with NO pipeline counterpart, named; whether anything still
   references `results/figures/`; and `reports/MANUSCRIPT_FIGURE_LINKAGE.md` written.
8. `results/`, `data/`, `manuscript/` (content), `config.R` and the submodule pointer byte-identical.
9. PR opened; `origin/main` unchanged.

CONSTRAINTS
- No analysis. No constant changed. No number changed. If making the pipeline run requires changing a
  number, STOP and report it as a finding.
- Reapply, do not cherry-pick. Adapt to his numbering rather than reverting it.
- Do not touch the manuscript beyond confirming it builds. PART D diagnoses and specifies; the fix
  is the next prompt. But it states a requirement, not a suggestion - do not soften it in the PR body.
- Do not touch `results/`, `data/`, or the `cauris_etcgem` submodule pointer.
- Keep the diff reviewable in ten minutes. Its value is that it is obviously safe.
- Autonomous; commit in parts: "packaging: pinned renv environment + 00_install + SETUP",
  "packaging: headless guards on the three Shiny scripts",
  "packaging: --file= resolution under the current numbering",
  "packaging: run_all.sh covering all eighteen stages, with timing and logging",
  "packaging: gitignore hygiene incl. LibreOffice lock files",
  "packaging: manuscript figure-linkage diagnosis and remediation spec".
```
