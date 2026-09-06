# Claude Code prompt — C13: recover the C10/C11 work stranded by an out-of-order merge (autonomous)

Run from the project root (`.../Candidas TPC/Candidas`). Recovery and reapplication. It changes no
analysis decision, no constant and no fitted number. It brings nine commits currently sitting on a
branch onto `main`, and reapplies one part of them to a layout that has since changed underneath.

NOTE TO USER: launch in an auto-approving mode. Moderate — a merge with a contained set of conflicts,
a mechanical reapplication, and a render. Tell Ilgaz before running it: it is his manuscript
flattening being preserved, and he may prefer to do the merge himself.

CONTEXT — WHAT HAPPENED AND WHY IT MATTERS.

Three Candidas pull requests were merged bottom-up but out of order. `#2` merged `gyd/packaging` into
`main` at 10:22:13; `#3` then merged into `gyd/packaging` at 10:22:54 and `#4` into
`gyd/figure-sourcing` at 10:23:39 — into parent branches that had already been merged away. The
result is that `origin/gyd/figure-sourcing` is 9 commits ahead of `main` and 4 behind, and none of
the C10 or C11 work reached `main`.

Verified on `origin/main` as it stands:

  * `scripts/config.R` contains no `fig_keep_add` (4 occurrences on the branch) — the whitelist fix
    is absent, so `18_n0_treatment_panel.R` still exits 0 and writes nothing;
  * `results/figures/Fig_n0_treatment_panel.png` does not exist — Supplementary Figure 6 still cannot
    be regenerated;
  * `reports/C11_scale_free/` does not exist;
  * the manuscript still embeds `media/media/image1..10.png`.

Meanwhile `b987f68` on `main` rewrote the manuscript to LEAD with the C11 scale-free result. So the
paper currently cites conclusions whose supporting analysis is not in the repository. That is the
failure mode C1 existed to catch, arrived at by a different route, and it is why this is urgent
rather than tidy-up.

THE MERGE IS FAVOURABLE. `main`'s four commits touch ONLY `manuscript/` and two PNGs under
`reports/review_figures/`. So of the nine stranded commits:

  * `117f861` (script 18 + `config.R` whitelist fix) — no conflict. `main` never touched `config.R`;
    its one change to script 18 is the `--file=` fix from our own packaging PR, in a different part
    of the file.
  * the four C11 commits (`reports/C11_scale_free/`) — no conflict, entirely new paths.
  * `c14241d`, `6d16262`, `7a46835` (manuscript figure sourcing) — CONFLICT. They edit
    `manuscript/draft/v3_quarto/*`, which `dbdb59d` relocated to a single canonical
    `manuscript/v3.qmd`.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: `git log --oneline
origin/main..origin/gyd/figure-sourcing` and the reverse; `reports/MANUSCRIPT_FIGURE_LINKAGE.md` on
the branch (the image-to-source mapping C9 produced); `manuscript/v3.qmd` and `_quarto.yml` on `main`
(the new flattened layout); and `dbdb59d` (what the flattening actually did). Understand both layouts
before merging anything.

PART A - merge, with one explicit resolution rule
- `git fetch origin`, then `git switch -c gyd/recover-figure-sourcing origin/main`.
- `git merge origin/gyd/figure-sourcing`.
- Everything outside `manuscript/` should merge clean. Report anything that does not — that would
  contradict the analysis above and needs understanding before proceeding.
- FOR EVERY CONFLICTED PATH UNDER `manuscript/`, TAKE MAIN'S SIDE. Git will present the branch's
  edits to `manuscript/draft/v3_quarto/*` as modify-versus-delete conflicts. Resolving them the wrong
  way RESURRECTS the directory Ilgaz deliberately retired in `dbdb59d`. Accept the deletions; accept
  main's `manuscript/v3.qmd`.
- After resolving, verify `manuscript/draft/` does not exist and that `manuscript/` matches `main`'s
  layout exactly. State that you have checked.
- Confirm the recovered content is present: `fig_keep_add` in `config.R`, `reports/C11_scale_free/`,
  and the script-18 fix.

PART B - confirm script 18 works on the merged tree
- Run `18_n0_treatment_panel.R` and confirm it now writes
  `results/figures/Fig_n0_treatment_panel.png`. It runs three nested Bayesian refits, so allow time.
- Confirm `09_bayesian_models.R` standalone output is unchanged by the whitelist fix — that was
  verified in C10 and should still hold after the merge; re-verify rather than assume.

PART C - reapply the figure sourcing to the new layout
- `manuscript/v3.qmd` currently embeds `media/media/image1..10.png`. Repoint each at its
  pipeline-produced file, using the mapping in `reports/MANUSCRIPT_FIGURE_LINKAGE.md`. C9 established
  that seven of nine were byte-identical to the pipeline output, so the substitution is mechanical
  rather than a rethink — but re-verify each against the CURRENT `results/figures/`, since the
  pipeline has been re-run since.
- Check the relative paths resolve from the new location. The old manuscript lived at
  `manuscript/draft/`; `v3.qmd` is now at `manuscript/`, so any path that worked before is now wrong
  by one directory level.
- Supplementary Figure 6 now has a producing script, so it sources like the rest. Note in the commit
  that the embedded version predated the current script — same values, but it carried per-panel
  "Clade IV = rank n" annotations the current script does not draw. Do not attempt to restore those
  annotations; that is Ilgaz's call and it is recorded in the C10 PR.
- Remove the embedded bitmaps ONLY after PART D verifies the render. If any image cannot be mapped to
  a pipeline source, leave it embedded and report it.

PART D - verify
- `quarto render` from wherever the flattened manuscript expects to be built. Report page count and
  figure count, and compare against the pre-change render: same pages, same figure order, nothing
  missing, no unresolved cross-reference.
- REDO THE PROPAGATION TEST. Regenerate one figure via its producing script, re-render, and confirm
  the change reaches the manuscript. Restore the original afterwards. This is the property the whole
  exercise exists to restore and it must be demonstrated, not assumed.
- Confirm `results/`, `data/`, `scripts/original_scripts/` and the `cauris_etcgem` submodule pointer
  are byte-identical to `origin/main`, by checksum, with the sole expected addition of
  `Fig_n0_treatment_panel.png`.

PART E - record what happened
- Add a short section to `reports/MANUSCRIPT_FIGURE_LINKAGE.md`, or a note beside it, recording the
  merge-order accident: which PRs merged when, what was stranded, and how it was recovered. Nobody
  will reconstruct this from the git history later, and it is the sort of thing that recurs.
- Do NOT edit any manuscript text. This prompt changes figure paths and nothing else.

PART F - PR
- Push and open a PR against `origin/main`. Title: "Recover C10/C11 work stranded by out-of-order
  merge; reapply figure sourcing to the flattened manuscript". Body: what was stranded and why, the
  resolution rule used for `manuscript/`, confirmation that script 18 now produces its figure, the
  figure-by-figure substitution, the propagation test result, and an explicit statement that no
  fitted number changes.
- Do not push to `main`. Do not force-push. Verify `origin/main` unchanged.
- Afterwards, `origin/gyd/figure-sourcing` and `origin/gyd/scale-free` are fully recovered and can be
  deleted — but confirm that by diffing against the merged result before recommending it, and leave
  the deletion to the user.

VERIFY (report all)
1. Which paths conflicted; confirmation all were under `manuscript/` and all resolved to main's side.
2. `manuscript/draft/` absent; `manuscript/` matches main's layout.
3. `fig_keep_add` present in `config.R`; `reports/C11_scale_free/` present; script-18 fix present.
4. Script 18 writes its figure; `09` standalone output unchanged, re-verified.
5. Figure-by-figure: which mapped byte-identically to the current `results/figures/`, which differed,
   which were left embedded and why.
6. Render: page and figure count against the pre-change render; anything missing or unresolved.
7. The propagation test result.
8. Checksums: `results/`, `data/`, submodule pointer unchanged bar the one expected addition.
9. PR opened; `origin/main` unchanged; whether the two branches are now fully recovered.

CONSTRAINTS
- No analysis decision, no constant, no fitted number, no manuscript TEXT edit.
- Take main's side for every conflict under `manuscript/`. Do not resurrect `manuscript/draft/`.
- If the merge conflicts anywhere outside `manuscript/`, STOP and report rather than resolving — that
  would mean the situation differs from what was analysed.
- Re-verify the figure mappings against the current `results/figures/` rather than trusting C9's
  record; the pipeline has been re-run since it was written.
- Autonomous; commit in parts: "merge: recover C10/C11 work stranded by out-of-order merge",
  "manuscript: source figures from results/ under the flattened layout",
  "manuscript: remove superseded embedded media",
  "docs: record the merge-order accident and its recovery".
```
