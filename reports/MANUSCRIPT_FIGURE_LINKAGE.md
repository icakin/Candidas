# Manuscript figure linkage — the record

**Status: RESOLVED.** Every figure in the manuscript resolves to a file produced by
the pipeline, referenced by path. The embedded bitmaps are gone. This document was
written as a diagnosis and remediation spec (C9); it is now the record of the final
state (C10), with the residual items kept at the end.

---

## 1. The requirement

**Every figure in the manuscript must resolve to a file produced by the pipeline,
referenced by path.** Embedded or extracted bitmaps are not an acceptable durable
arrangement, whatever the short-term convenience of a Word round-trip, because:

* **A regenerated figure never reaches the paper.** The manuscript can then silently
  disagree with the analysis it reports — the failure this project has already hit
  once, and the reason the C-series exists.
* **A referee cannot trace a figure back to the code that made it.** Traceability from
  claim → figure → script → data is the property the whole exercise buys.

## 2. Final state

`manuscript/v3.qmd` references six figures, all by path into `results/`.
`manuscript/draft/v3_quarto/media/` has been deleted.

| # | Manuscript figure | Path referenced from `v3.qmd` | Producing script |
|---|---|---|---|
| 1 | Figure 1 — growth/respiration decoupling | `../results/figures/manuscript/FIG1_decoupling.png` | `16_fig1.R / 17_fig2.R` |
| 2 | Figure 2 — the carbon cost of fever | `../results/figures/manuscript/FIG2_consequences.png` | `16_fig1.R / 17_fig2.R` |
| 3 | Supp. Fig. 1 — equilibration sensitivity | `../results/figures/Fig_temperature_equilibration.png` | `08_temperature_equilibration_sensitivity.R` |
| 4 | Supp. Fig. 2 — respiration + envelope | `../results/figures/fig_bayes_resp_arrhenius.png` | `11_bayesian_plots.R` |
| 5 | Supp. Fig. 3 — CUE + envelope | `../results/figures/fig_bayes_cue_by_clade.png` | `11_bayesian_plots.R` |
| 6 | Supp. Fig. 4 — N0 treatment panel | `../results/figures/Fig_n0_treatment_panel.png` | `15_n0_treatment_panel.R` |

**Built but not yet placed.** Two finished main-figure candidates exist as files with
their own generators and captions, and are deliberately NOT referenced by `v3.qmd` while
their placement is undecided:

| Figure | File | Generator | Caption |
|---|---|---|---|
| thermal/phylogeny discordance | `results/figures/manuscript/FIG3_discordance.png` | `scripts/18_fig3.R` | `gem/FIG3_v8_caption.md` |
| etcGEM counterfactual | `results/figures/manuscript/FIG4_etcgem_counterfactual.png` | `scripts/19_fig4.R` | `gem/FIG4_etcgem_caption.md` |

A draft of `v3.qmd` carrying both, with their captions and placeholder Results text, is
kept at `manuscript/v3_withfig34_20260904.qmd`. Adopting it is a single file copy; until
then this table describes what the manuscript actually references, which is the point of
the document.

**Superseded numbering.** Earlier revisions of this document listed the etc-GEM pipeline
schematic (`FIG_model_schematic.png`, `17_schematic.py`) as Figure 3 and the capacity axis
(`FIG_MODEL.png`, `16_fig1.R / 17_fig2.R`) as Figure 4, and listed a ten-figure manuscript.
Neither was in fact referenced by `v3.qmd`; the document had drifted from the file it
describes. Those two figures remain unreferenced and should be placed in the supplement or
dropped deliberately rather than left in this table.

**A broken reference this table did not catch.** Until 2026-09-04, `v3.qmd` referenced
`FIG2_the_bill.png`, which has never existed in the repository. `16_fig1.R / 17_fig2.R` writes
`FIG2_consequences` via `save_fig()`, but its own OUTPUTS header comment named
`FIG2_the_bill.png`, and the manuscript followed the comment. Figure 2 therefore could not
render. Both the manuscript path and the stale comment are corrected, and the comment now
says explicitly that the name follows `save_fig()` and not the header. The lesson for this
document is that listing a path is not the same as checking it resolves; the check is one
line and is now recorded in §6.

**Two directories, deliberately.** `14`, `15` and `16` write into
`results/figures/manuscript/`; `08`, `11` and `18` write into `results/figures/`.
Copying the latter into `manuscript/` would tidy the paths at the cost of a second
copy of each figure — and the existing top-level/`manuscript/` duplicates have already
drifted apart (§5). One figure, one path.

### 2.1 Propagation is demonstrated, not assumed

Regenerating a figure through its own producing script and re-rendering:

```
BEFORE regeneration: (2082, 2762)  sha=0170848321a57bde
$ python3 archive/old_etcgem/17_schematic.py
AFTER  regeneration: (2080, 2762)  sha=3ab458b8f3095bac
$ cd manuscript/draft && quarto render
page 6 of the rendered PDF: 2080x2762
```

The change reached the paper. Before this, it could not have.

## 3. Supplementary Figure 6: what was wrong and what fixed it

`15_n0_treatment_panel.R` exited 0, printed `wrote Fig_n0_treatment_panel.png
(Supplementary Figure 6)`, and wrote nothing at all. It was the one bitmap of ten
with no pipeline counterpart, because its only producing script could not produce it.

**The mechanism.** `config.R` gates `ggsave` on a whitelist. `18` appended its figure
name to the `FIG_KEEP` global to get past the gate, then called `09_bayesian_models.R`
three times; `09`'s `source(file.path(.this_dir, "config.R"))` defaults to
`local = FALSE`, so it re-evaluates **in the global environment regardless of who
called it**, resetting `FIG_KEEP` to its base list. Reproduced directly:

```
after 18 line 59        : FIG_KEEP has 30 names, gate allows the panel = TRUE
after one nested 09     : FIG_KEEP has 29 names, gate allows the panel = FALSE
ggsave() returned       : NULL (write suppressed)
file written            : FALSE
```

Suppressed at `config.R:533`, inside `.gated_ggsave`.

**The fix.** The extension point is no longer a clobberable global. `config.R` gains
`fig_keep_add()`, which records extra names in an **option** that `config.R` never
resets; `.fig_keep_ok()` tests `FIG_KEEP` union that option. `18` now calls
`fig_keep_add("Fig_n0_treatment_panel")`. Any future script that invokes `09` gets the
same guarantee without having to know the hazard exists. `.gated_ggsave` and
`.gated_pdf` now `message()` when they drop a write — the silence is what hid this.

*Rejected:* making `09` source `config.R` with `local = TRUE`. It is the more general
root cause, but it relocates ~100 config objects for three callers and risks changing
a number. Kept as a residual item (§6).

**`09` is unchanged by the fix.** Run standalone before and after, all 14 outputs are
byte-identical.

**The regenerated figure shows the same thing as the embedded one.** The table behind
it, `results/tables/bayes_resp_arr_E_three_treatments.csv`, came back **byte-identical**
to the committed copy: max |E difference| = **0.00000 eV** across all 15
clade × treatment values, rank vectors identical (Clade IV = 1, 1, 4).

Rendering differences, all cosmetic and all reported rather than glossed:

| | embedded `image10.png` | regenerated |
|---|---|---|
| size | 3560 × 1004 | 2400 × 720 (`ggsave` 12 × 3.6 in at 200 dpi) |
| title / subtitle | absent | present |
| per-panel "Clade IV = rank *n*" annotation | present | **absent** |
| x-axis range | to 0.8 | to ~0.7 |

The embedded bitmap therefore came from an **earlier revision of `18`** whose plotting
code has since changed. The values, clades and treatments are identical, so this is not
material — and the dropped annotation states a fact the caption already states in prose
("Clade IV holds the lowest activation energy under the published and
equilibration-corrected treatments and moves to fourth of five when the term is
removed"), which the regenerated ranks confirm. Nothing in the caption refers to an
in-figure annotation. This drift is exactly what sourcing by path prevents in future.

## 4. Render verification

| | before | after |
|---|---|---|
| pages | 14 | 14 |
| figures | 10 | 10 |
| figure pages | 2, 4, 6, 8, 11, 12, 12, 13, 13, 14 | identical |
| missing figures / unresolved refs | — | none |

12 of the 14 pages are **pixel-identical** to the pre-change render. The three that
move, and why:

| page | mean abs. pixel diff | why |
|---|---|---|
| 2 (Figure 1) | 0.218 | now the live `manuscript/` copy instead of the stale top-level one; below this figure's own 0.538 run-to-run jitter |
| 12 (Supp. Fig. 3) | 0.747 | lossy Word re-encode replaced by the authoritative file |
| 14 (Supp. Fig. 6) | 3.220 | by design — the figure is now produced by the pipeline |

## 5. The duplicate figures under `results/figures/`

`16_fig1.R / 17_fig2.R`'s `save_fig()` writes **only** to `results/figures/manuscript/`.
The top-level copies of `FIG1_decoupling.png`, `FIG2_the_bill.png`, `FIG_MODEL.png` and
`FIG_MODEL_SUPP_validation.png` are therefore **not** rewritten by the pipeline — they
are older files, and they differ from the `manuscript/` copies (mean abs. pixel
difference 0.45, 0.38, 1.20, 1.76 respectively). `FIG_model_schematic.png` and
`FIG_MODEL_SUPP_consistency.png` are identical in both places.

The manuscript now references the live `manuscript/` copies throughout. The stale
top-level duplicates should be deleted so there is exactly one file per figure — not
done here, because this change touches no committed result.

## 6. Residual items

None of these blocks anything; all are cheap and worth doing.

1. **`09` still sources `config.R` into the global environment.** The figure whitelist
   is now immune, but any other caller state named like a config variable is still
   clobbered by a nested `09`. The general fix is `local = TRUE`, which needs its own
   before/after byte-identity check across `09`, `10`, `18` and `run_all.R`.
2. **The figures are not reproducible run to run.** Two identical `run_all.sh` runs
   give byte-identical tables but different PNGs — `boxplot_CUE.png` mean 1.21,
   `boxplot_K_O2_rate.png` 1.00, `manuscript/FIG1_decoupling.png` 0.34. Unseeded
   jitter in the plotting layer, plus `12_carbon_tax.R`'s unseeded
   `sample.int(nrow(p), 1500)`, which also makes `results/tables/carbon_tax_curves.csv`
   the one table that moves between identical runs. Seed both, and "the figure changed"
   starts to mean something.
3. **`17_schematic.py` uses `bbox_inches="tight"`,** so the schematic's canvas size is
   environment-dependent (2080 × 2762 here vs 2082 × 2762 committed). Pin an explicit
   `figsize` if it is to be checksummable.
4. **`v3.tex` is a tracked build artefact** that every render rewrites. It is refreshed
   in this change so it no longer points at deleted files, but it will drift again.
   Add it to a `manuscript/.gitignore` and `git rm --cached` it. The same now
   applies to the tracked `manuscript/v3.pdf`.
5. **Stale script-number references in comments,** left over from the renumbering:
   `10_n0_term_test.R:95` says "re-run 08" where it means 09;
   `13_uncertainty_bands.R:7` and `:17` say "from script 11" where they mean 08.
6. **A regression guard.** A check that greps `v3.qmd` for `media/media/` and fails if
   anything matches, plus an assertion that every path it references exists under
   `results/`. Without something of this shape the arrangement decays the next time the
   manuscript goes through Word.


---

## 7. The merge-order accident, and how the work was recovered

Recorded because nobody will reconstruct this from the git history later, and
because it is the sort of thing that recurs.

### What happened

Three pull requests were stacked — packaging → figure sourcing → scale-free
analysis — and merged bottom-up, but **out of order**:

| PR | merged | into |
|---|---|---|
| #2 `gyd/packaging` | 10:22:13 | `main` |
| #3 `gyd/figure-sourcing` | 10:22:54 | `gyd/packaging` |
| #4 `gyd/scale-free` | 10:23:39 | `gyd/figure-sourcing` |

`#2` landed on `main` **41 seconds before** `#3` merged into the branch `#2` had
just consumed. From that moment `gyd/packaging` and `gyd/figure-sourcing` were
parent branches that had already been merged away, and everything subsequently
merged into them accumulated off to one side. `origin/gyd/figure-sourcing`
finished 9 commits ahead of `main` and 4 behind, and none of the C10 or C11 work
ever reached `main`.

### What was stranded

* `scripts/config.R` had no `fig_keep_add`, so `15_n0_treatment_panel.R` still
  exited 0 and wrote nothing — Supplementary Figure 6 could not be regenerated;
* `results/figures/Fig_n0_treatment_panel.png` did not exist;
* `reports/C11_scale_free/` did not exist;
* the manuscript still embedded `media/media/image1..10.png`.

**Why that was urgent rather than tidy-up.** Commit `b987f68` on `main` had
already rewritten the manuscript to *lead* with the C11 scale-free result. So the
paper cited conclusions whose supporting analysis was not in the repository —
the same class of failure the C-series exists to catch, arrived at from the
opposite direction. The analysis was right and the paper was right; the audit
trail between them was missing.

### How it was recovered

`gyd/figure-sourcing` was merged into a branch off `main`. The merge was
favourable because `main`'s four commits touch only `manuscript/` and two renames
under `reports/review_figures/`:

* **ten conflicts, all rename/delete under `manuscript/media/media/`**, and none
  anywhere else. Every one resolved to **main's** side, preserving the flattened
  layout of `dbdb59d`. Resolving them the other way would have resurrected the
  `manuscript/draft/` nesting that commit deliberately retired.
* `manuscript/v3.qmd` and `v3.tex` **auto-merged** rather than conflicting, and
  the auto-merge silently took the branch's `../../../results/...` paths while
  keeping main's text. Those were correct from `manuscript/draft/v3_quarto/` and
  are wrong by two levels from `manuscript/`. Both were reset to main's version
  in the merge commit and repointed properly afterwards, as `../results/...`.

That second point is the trap worth remembering: **a clean auto-merge is not
evidence of a correct one.** Relative paths merge textually and break silently.

### The lesson

Merge a stack **top-down**, or rebase each PR onto the new base after its parent
lands. Merging bottom-up works only if every merge happens before its parent is
consumed, which is a race nobody should have to win. If GitHub's "merge" button
is used on a stack, check afterwards that the branch you merged *into* is still
an ancestor of `main`.


---

## 6. Path check

Every path referenced by the manuscript must resolve. This is the check, and it is the one
that would have caught the Figure 2 break on the day it was introduced:

```
cd manuscript && grep -o '(\.\./results/figures[^)]*)' v3.qmd | tr -d '()' \
  | while read p; do [ -f "$p" ] && echo "OK      $p" || echo "MISSING $p"; done
```

As of 2026-09-04 all six resolve.
