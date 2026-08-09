# Manuscript figure linkage — diagnosis and remediation spec

**Status:** diagnosis complete, remediation NOT applied. This document exists so the
fix is a small specified job rather than an investigation. No manuscript file is
changed by the pull request that adds this document.

---

## 1. The requirement

**Every figure in the manuscript must resolve to a file produced by the pipeline,
referenced by path.** Embedded or extracted bitmaps are not an acceptable durable
arrangement, whatever the short-term convenience of a Word round-trip.

This is a requirement, not a preference, because two things break when the link is
severed:

* **A regenerated figure never reaches the paper.** The manuscript can then silently
  disagree with the analysis it reports. That is not hypothetical here: it is the
  failure this project has already hit once, and the reason the C-series exists.
* **A referee cannot trace a figure back to the code that made it.** Traceability
  from claim → figure → script → data is the property the whole exercise buys. A
  bitmap with no path breaks the chain at its first link.

## 2. Current state

`manuscript/draft/` was rebuilt for v3. The hand-written `_body_main.qmd`,
`_body_supp.qmd`, `manuscript.qmd` and `supplementary.qmd` are gone, replaced by
`manuscript/draft/v3_quarto/` containing `v3.qmd`, `v3.tex` and ten **extracted
bitmaps** at `media/media/image1..10.png`.

**Nothing under `manuscript/` references `results/figures/` any more.** Verified with:

```bash
grep -rn "results/" manuscript/ --include='*.qmd' --include='*.tex' \
                                --include='*.yml' --include='*.md'
# (no matches)
```

Every one of the ten figure references is of the form
`![](media/media/imageN.png){width="..." height="..."}`. There is **no live pipeline
figure input to the manuscript**: the build reads only the extracted bitmaps.

## 3. Inventory: every embedded image → its pipeline figure → its producing script

Matched by pixel dimensions and by whole-file SHA-256 against the committed
`results/figures/` tree (see `§3.1` for method).

| Embedded | px | Manuscript figure | Pipeline file | Match | Producing script |
|---|---|---|---|---|---|
| `image1.png` | 4322×5078 | Figure 1 — growth/respiration decoupling | `results/figures/FIG1_decoupling.png` | **byte-identical** | `14_main_figures.R` (`save_fig(FIG1, "FIG1_decoupling", 215)`) |
| `image2.png` | 4322×4960 | Figure 2 — the carbon cost of fever | `results/figures/manuscript/FIG2_the_bill.png` | **byte-identical** | `14_main_figures.R` (`save_fig(FIG2, "FIG2_the_bill", 210)`) |
| `image3.png` | 2082×2762 | Figure 3 — etc-GEM pipeline schematic | `results/figures/FIG_model_schematic.png` (and the `manuscript/` copy — the two are identical) | **byte-identical** | `17_schematic.py` |
| `image4.png` | 4322×5078 | Figure 4 — capacity axis | `results/figures/manuscript/FIG_MODEL.png` | **byte-identical** | `14_main_figures.R`, Fig-3/etc-GEM section (`save_fig(FIG, "FIG_MODEL", 215)`) |
| `image5.png` | 4322×3543 | Supplementary Figure 1 — capacity vs coding variation | `results/figures/manuscript/FIG_MODEL_SUPP_validation.png` | **byte-identical** | `16_supplementary_figures.R` (line ~331) |
| `image6.png` | 2125×1653 | Supplementary Figure 2 — self-consistency | `results/figures/FIG_MODEL_SUPP_consistency.png` (and the `manuscript/` copy — identical) | **byte-identical** | `14_main_figures.R` (`save_fig(fs, "FIG_MODEL_SUPP_consistency", 70)`) |
| `image7.png` | 1950×690 | Supplementary Figure 3 — equilibration sensitivity | `results/figures/Fig_temperature_equilibration.png` | same size and content, **not byte-identical** (mean abs. pixel difference 6.06/255) | `08_temperature_equilibration_sensitivity.R` |
| `image8.png` | 3900×2250 | Supplementary Figure 4 — respiration + equilibration envelope | `results/figures/fig_bayes_resp_arrhenius.png` | **byte-identical** | `11_bayesian_plots.R` (line ~225) |
| `image9.png` | 3900×2250 | Supplementary Figure 5 — CUE + equilibration envelope | `results/figures/fig_bayes_cue_by_clade.png` | **byte-identical** | `11_bayesian_plots.R` (line ~270) |
| `image10.png` | 3560×1004 | **Supplementary Figure 6 — N0 treatment panel** | **NONE** | **no pipeline counterpart** | see §4 |

### 3.1 Method

```
for each image1..10.png:
    exact  : sha256(embedded) == sha256(candidate) over results/figures/**/*.png
    content: same (w,h) -> mean and max |pixel difference| in RGB
```

`image7` is the only matched figure that is not byte-identical: it is the same plot
(same titles, same three panels, same curves) re-encoded by the Word round-trip. That
is itself evidence the round-trip is lossy — a figure that survives it *looking* right
is not the same file, and cannot be checksummed against the pipeline.

### 3.2 A complication the remediation has to know about: the figures are not reproducible

Two *identical* back-to-back runs of `bash scripts/run_all.sh` on the same machine
produce byte-identical **tables** (one exception, below) but **different PNGs**:

| figure | run 2 vs run 1, mean abs. pixel difference |
|---|---|
| `boxplot_CUE.png` | 1.21 |
| `boxplot_K_O2_rate.png` | 1.00 |
| `manuscript/FIG1_decoupling.png` | 0.34 |
| `manuscript/FIG2_the_bill.png` | 0.11 |
| `manuscript/FIG_MODEL_SUPP_validation.png` | 0.15 |
| `FIG_carbon_tax.png` | 0.16 |

These are the same order as, and in several cases larger than, the difference between
a regenerated figure and its committed copy. The cause is unseeded RNG in the plotting
layer (jittered point positions), plus `12_carbon_tax.R`'s unseeded
`sample.int(nrow(p), 1500)`, which also makes `results/tables/carbon_tax_curves.csv`
the one table that does not reproduce run to run.

This matters for the remediation in two ways:

* **It strengthens the requirement, it does not weaken it.** With bitmaps embedded,
  there is no way to tell a real change from jitter — you cannot checksum the paper
  against the pipeline at all. With paths, at least the file identity is exact.
* **Seed the plots first.** Before wiring the `.qmd` to `results/figures/`, set a seed
  for the jitter and for `12`'s subsample, so "the figure changed" means something.
  Until then, expect a rendered page to differ cosmetically on every run.

One more mismatch of the same family: `17_schematic.py` regenerates
`FIG_model_schematic.png` at **2080×2762**, while the committed copy is **2082×2762**.
matplotlib's `bbox_inches="tight"` depends on font metrics, so the schematic's canvas
size is environment-dependent. Pin it with an explicit `figsize`/`bbox_inches` rather
than `tight` if the figure is to be checksummable.

## 4. The serious case: `image10.png` (Supplementary Figure 6)

**`image10.png` has no pipeline counterpart. There is no file under `results/` that it
matches, at any size.**

What it shows: three side-by-side panels of the hierarchical Bayesian respiration
activation energy `E_R` per clade (posterior median, 95% CrI), one panel per treatment
of the inoculum back-projection — *"With term (published)"*, *"Equilibration-corrected"*
and *"Term-free"* — each annotated with Clade IV's rank (1, 1 and 4 respectively).
Clade IV is drawn in red.

That is exactly the figure `scripts/18_n0_treatment_panel.R` is written to produce
(`Fig_n0_treatment_panel.png`, its own header calls it "Supplementary Figure 6"). So a
producing script does exist — but:

1. **Its output is not committed.** `results/figures/Fig_n0_treatment_panel.png` does
   not exist in the repository. Every other manuscript figure has a committed copy.
2. **`18` is not run by any runner.** It is not in `run_all.R`, and `run_all.sh`
   deliberately leaves it out: it refits the Bayesian respiration model up to three
   times and overwrites `results/tables/bayes_resp_arr_summary.csv` and `results/rds/`
   with alternative fits without restoring them.
3. **The committed bitmap does not match what the script would now write.** `18` calls
   `ggsave(..., width = 12, height = 3.6, dpi = 200)` → 2400×720 px, and its `ggplot`
   carries a `title` and a `subtitle`. `image10.png` is 3560×1004 px and has neither.
   The embedded figure is therefore from a *different* version of this plot than the
   one in the tree.

### 4.1 It is worse than "not committed": the script cannot write the figure

`18_n0_treatment_panel.R` was run to completion under the pinned environment. It
**exits 0 and prints `wrote Fig_n0_treatment_panel.png (Supplementary Figure 6)`** —
and the file does not exist, anywhere on disk. The write is silently dropped.

The cause is `config.R`'s figure whitelist, and `18` knows about it. Line 52:

```r
if (exists("FIG_KEEP")) FIG_KEEP <- unique(c(FIG_KEEP, "Fig_n0_treatment_panel"))
```

That works — for as long as it lasts. `18` then sources `09_bayesian_models.R` three
times, and `09` begins with `source(file.path(.this_dir, "config.R"))`. `source()`
defaults to `local = FALSE`, which evaluates **in the global environment regardless of
who called it**, so each nested `09` re-runs `config.R`'s top-level
`FIG_KEEP <- c(...)` and resets the list to its 29 original names. By the time `18`
reaches its `ggsave()`, the gate rejects the filename. Demonstrated directly:

```
AFTER 18's append            : in FIG_KEEP = TRUE  | gate says write = TRUE
AFTER one nested 09 -> config: in FIG_KEEP = FALSE | gate says write = FALSE
```

The gated `ggsave` returns `invisible(NULL)` on rejection, and `18`'s success message
is unconditional, so nothing reports the failure.

So Supplementary Figure 6 is, today, a claim in the paper whose only producing script
**cannot produce it**, and says it did. It is the same failure mode this project has
already hit once, and it is the item on this list that needs fixing first.

Two ways to fix it, in `18`:

* re-append to `FIG_KEEP` *after* the last nested `09` call and immediately before
  `ggsave()` — smallest change, still fragile; or
* better, have the gate consult an option rather than a global that any `source()` can
  clobber, e.g. `options(candidas.fig_keep_extra = ...)` read inside `.fig_keep_ok()`.

Whichever is chosen, `.gated_ggsave` should `message()` when it drops a write. A
silent no-op is how this got here.

## 5. Pipeline figures that are NOT in the manuscript

Not a defect — recorded so the remediation does not accidentally wire them in:

* `results/figures/FIG_MODEL_SUPP.png` (4322×4488, from `16_supplementary_figures.R`)
* `results/figures/FIG_carbon_tax.png` (from `12_carbon_tax.R`)
* `results/figures/fig_cue_uncertainty.png`, `fig_respiration_uncertainty.png`
  (1500×975, from `15_uncertainty_bands.R`) — superseded for the manuscript by the
  envelopes `11_bayesian_plots.R` now draws directly onto Supp. Figs 4 and 5
* `fig_bayes_contrasts.png`, `fig_bayes_cue_by_isolate.png`,
  `fig_bayes_growth_tpc_by_clade.png`, `fig_bayes_E_resp_minus_growth.png`, the
  `boxplot_*` set, `CUE_vs_T.png`, `arrhenius_respiration_fgC_h.png`,
  `TPC_*_by_isolate.pdf`

## 6. A second, smaller inconsistency to fix at the same time

`14_main_figures.R`'s `save_fig()` writes **only** to `results/figures/manuscript/`. The
top-level copies `results/figures/FIG1_decoupling.png`, `FIG2_the_bill.png`,
`FIG_MODEL.png` and `FIG_MODEL_SUPP_validation.png` are therefore **not** rewritten by
the current pipeline — they are older files left in place, and they differ from the
`manuscript/` copies:

| File | top-level vs `manuscript/` copy |
|---|---|
| `FIG1_decoupling.png` | differ (mean abs. pixel difference 0.45) |
| `FIG2_the_bill.png` | differ (0.38) |
| `FIG_MODEL.png` | differ (1.20) |
| `FIG_MODEL_SUPP_validation.png` | differ (1.76) |
| `FIG_model_schematic.png` | identical |
| `FIG_MODEL_SUPP_consistency.png` | identical |

The manuscript currently mixes them: `image1` came from the **top-level** copy while
`image2`, `image4` and `image5` came from the **`manuscript/`** copy. Whichever
directory the `.qmd` is pointed at, the stale duplicates in the other should go, so
there is exactly one file per figure and it is the one the pipeline writes.

## 7. Remediation spec

Everything below is mechanical. It is a follow-up prompt, not a research task.

### 7.1 Give Supplementary Figure 6 a producing, committed output

1. Decide what Supp. Fig. 6 should look like: `18_n0_treatment_panel.R` currently adds
   a title and subtitle that the embedded bitmap does not have. Either drop them from
   the script (the caption in `v3.qmd` already carries that text) or accept them.
2. Run `Rscript scripts/18_n0_treatment_panel.R`, then re-run
   `09 → 11 → 12 → 14 → 15` to restore the published `bayes_*` tables and `results/rds/`
   that `18` overwrites.
3. Commit `results/figures/manuscript/Fig_n0_treatment_panel.png`.
4. Better: make `18` write into `results/figures/manuscript/` via the same `save_fig()`
   contract as the other manuscript figures, and make it non-destructive (fit into a
   scratch table rather than over `bayes_resp_arr_summary.csv`) so it can join
   `run_all.sh` as a real stage instead of a hand-run diagnostic. Until it is
   non-destructive it must stay out of the unattended runner.

### 7.2 Point the `.qmd` at pipeline paths

The build's working directory is `manuscript/draft/` (that is where `_quarto.yml`
lives and where `quarto render` is invoked). `v3.qmd` sits one level deeper, in
`manuscript/draft/v3_quarto/`. Quarto resolves image paths relative to the `.qmd`, so
from `v3.qmd` the results tree is `../../../results/figures/manuscript/`.

Replace, in `manuscript/draft/v3_quarto/v3.qmd`:

| Line | from | to |
|---|---|---|
| 36 | `media/media/image1.png` | `../../../results/figures/manuscript/FIG1_decoupling.png` |
| 48 | `media/media/image2.png` | `../../../results/figures/manuscript/FIG2_the_bill.png` |
| 58 | `media/media/image3.png` | `../../../results/figures/manuscript/FIG_model_schematic.png` |
| 66 | `media/media/image4.png` | `../../../results/figures/manuscript/FIG_MODEL.png` |
| 100 | `media/media/image5.png` | `../../../results/figures/manuscript/FIG_MODEL_SUPP_validation.png` |
| 104 | `media/media/image6.png` | `../../../results/figures/manuscript/FIG_MODEL_SUPP_consistency.png` |
| 108 | `media/media/image7.png` | `../../../results/figures/Fig_temperature_equilibration.png` |
| 112 | `media/media/image8.png` | `../../../results/figures/fig_bayes_resp_arrhenius.png` |
| 116 | `media/media/image9.png` | `../../../results/figures/fig_bayes_cue_by_clade.png` |
| 120 | `media/media/image10.png` | `../../../results/figures/manuscript/Fig_n0_treatment_panel.png` (after §7.1) |

Three of these (`image7`, `image8`, `image9`) live at the top level of
`results/figures/` rather than in `manuscript/`, because `08` and `11` write there.
Either accept the two directories, or add a `save_fig()`-style manuscript copy to `08`
and `11`. Prefer the latter, so the `.qmd` has exactly one figure directory.

Keep the existing `{width= height=}` attributes: they are the sizes the layout was
tuned for and are independent of where the file comes from.

**Caveat to check on the first render.** These paths climb out of the Quarto project
directory (`manuscript/draft/`). Quarto resolves them for LaTeX and copies them for
HTML/docx, but it has historically been fussy about resources above the project root.
If `quarto render` complains, the fix is one of: add
`resources: ["../../../results/figures"]` to `_quarto.yml`; move the project root up
to the repository root; or have `08`, `11`, `12`, `14`, `16` and `17` all write a copy
into a single `manuscript/draft/figures/` directory that the pipeline owns. The last
option is the most robust and is the one to reach for if the first render fights back
— it keeps a real path, produced by the pipeline, and puts it inside the project.

### 7.2a `v3.tex` is a tracked build artefact, and it is already stale

`manuscript/draft/v3_quarto/v3.tex` is generated by `keep-tex: true`, but it is
committed — so **every `quarto render` modifies a tracked file**. Worse, the committed
copy is out of date with the committed `.qmd`: it still attributes Supplementary
Figure 6 to `scripts/16b_n0_treatment_panel.R`, the name that script had before the
renumbering. Re-rendering corrects it to `scripts/18_n0_treatment_panel.R`.

Add `v3.tex` to `manuscript/draft/.gitignore` and `git rm --cached` it, so a render
stops dirtying the working tree and nobody reads a stale `.tex` as the source of
truth. (Not done in the packaging PR: it is a manuscript file, and that PR touches
none.)

### 7.3 Verify the render is unchanged

This is the check that makes the change safe, and it must be part of the same commit's
verification:

```bash
cd manuscript/draft && quarto render                 # before, on the current tree
cp -R _output /tmp/render_before
# apply 7.1 + 7.2
cd manuscript/draft && quarto render
# compare page images, not PDF bytes (timestamps and object ids differ):
for f in /tmp/render_before/*.pdf; do
  b=$(basename "$f")
  magick -density 100 "$f" /tmp/before_%03d.png
  magick -density 100 "_output/$b" /tmp/after_%03d.png
  # then compare each page pair; expect an empty difference
done
```

Expect: every page identical except Supplementary Figure 6, which changes by design
(§7.1) and must be checked by eye against the caption.

### 7.4 Then delete the bitmaps

Once §7.3 passes, `git rm -r manuscript/draft/v3_quarto/media/`. Leaving them in place
guarantees the next round-trip re-severs the link.

## 8. Guard against regression

Add a check that fails loudly rather than silently:

* a test that greps `v3.qmd` for `media/media/` and fails if any match survives, and
* a `make`-style step (or a stage in `run_all.sh`) that asserts every path referenced
  by `v3.qmd` exists under `results/` **and** is newer than the derived table it
  depends on.

Without something of this shape, the arrangement decays the next time the manuscript
goes through Word.
