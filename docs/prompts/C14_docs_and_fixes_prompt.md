# Claude Code prompt — C14: bring the C1/C2 reports and the prompt series onto main, plus two small fixes (autonomous)

Run from the project root (`.../Candidas TPC/Candidas`). Documentation, plus two contained fixes. It
runs no analysis, changes no constant and alters no fitted number. It closes an omission: the C1 and
C2 reports and the entire prompt series were excluded from the packaging PR to keep it reviewable,
and the follow-up that should have brought them across was never written.

NOTE TO USER: launch in an auto-approving mode. Small — file recovery, two edits, one PR. Main is
healthy now (`fig_keep_add` present, `reports/C11_scale_free/` present), so this is additive.

CONTEXT. `C9` deliberately excluded `reports/`, `prompts/`, `runs/` and `HANDOVER.md` so the packaging
PR stayed small. The OxygenModel repository got its equivalent material back through a separate docs
PR; Candidas never did. Confirmed absent from `origin/main` today:

    reports/REPRODUCTION.md      the C1 reproduction report
    reports/N0_SENSITIVITY.md    the C2 three-arm N0 report
    reports/n0/                  its rendered PDF and .qmd
    reports/RUNBOOK.md           from-scratch instructions
    prompts/                     the whole C-series specification record
    HANDOVER.md

All of it is on `origin/gyd/c1-c2-review` (PR #1, still open). `main` currently has only
`reports/C11_scale_free/`, `reports/MANUSCRIPT_FIGURE_LINKAGE.md` and `reports/review_figures/`.

This matters beyond tidiness: the manuscript's Robustness section rests on the C2 N0 work, and that
report is not in the repository.

TWO SMALL FIXES RIDE ALONG, both reported from the field:
  * `00_install.R`'s compile scan flags `lattice`, `Matrix` and `survival` as missing. They are
    base-recommended packages living in the SYSTEM library, so the scan is searching the wrong tree.
  * `renv.lock` records the Bioconductor repository as an `r-universe.dev` URL that serves none of the
    three Bioc packages as arm64 binaries, while `bioconductor.org/packages/3.22/bioc` has all three
    at exactly the pinned versions. This has been approved in principle, but the diff must be shown
    before it merges — so make the change, and present it clearly in the PR body.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: on
`origin/gyd/c1-c2-review` — `reports/`, `prompts/`, `HANDOVER.md`; on `origin/main` — the current
`reports/`, `scripts/00_install.R`, `renv.lock` (the Repositories block and the three Bioconductor
package entries), and `README.md`. Branch: `git switch -c gyd/docs origin/main`.

PART A - bring the documentation across
- Take `reports/REPRODUCTION.md`, `reports/N0_SENSITIVITY.md`, `reports/n0/`, `reports/RUNBOOK.md`,
  `prompts/` and `HANDOVER.md` from `origin/gyd/c1-c2-review`, plus the prompts written since that
  branch was cut (C9-C13 and this one) if they are not already tracked.
- Do NOT bring `runs/`. Those trees are large and non-canonical. If the C1/C2 reports read from
  `runs/`, apply the policy the OxygenModel repository settled on: track the small TEXT TABLES the
  reports actually read, exclude figures and `.rds`, and record the split in a manifest. Establish
  what they read by parsing the sources, not by guessing from directory names.
- Verify `git diff --stat origin/main...gyd/docs` shows only `reports/`, `prompts/`, `HANDOVER.md`,
  `.gitignore`, `README.md`, `scripts/00_install.R` and `renv.lock` — nothing else.

PART B - make the recovered documents true on today's main
They were written before the script renumbering, the manuscript flattening and the etc-GEM cut.
Correct stale DOCUMENTATION; annotate, never rewrite, superseded FINDINGS.
- `prompts/README.md`: mark C1, C2, C2b, C9-C13 as run with a one-line outcome each; note that C3-C8
  were described but never written, and that the etc-GEM prompts (C5-C8) are now moot because that
  section has been cut from the manuscript.
- `reports/RUNBOOK.md`: update for the CURRENT eighteen-script numbering, `run_all.sh`, and the
  flattened `manuscript/`. The old runbook predates all three.
- `HANDOVER.md`: refresh the status section, or retire it with a note pointing at the reports if
  `reports/README.md` now serves that purpose better. State which you chose and why.
- Add a `reports/README.md` if absent: one line per report — ID, title, date, conclusion.
- ADD DATED NOTES where a conclusion has been superseded. At minimum:
  * `REPRODUCTION.md` reported that the etc-GEM outputs did not reproduce (capacity regenerating
    1.33-1.37x higher, an anchor change, a temperature grid the code could not emit). That section is
    now historical: the etc-GEM has been CUT from the manuscript. Note that, and that the code and
    outputs remain in the repository with `scripts/README.md` recording what a future treatment would
    need.
  * `N0_SENSITIVITY.md` predates the in-code `Ninoc` regeneration and the C11 scale-free result. Note
    that the reconciliation subsequently closed (FC_Initial against N0 at fit start: 1.53x, range
    0.79-3.24) and that the CUE optimum is now established scale-free.
  Each report stands as the record of what was found at the time. The note stops a later reader
  taking its numbers as current.
- Add a short section recording the merge-order accident and its recovery, if
  `MANUSCRIPT_FIGURE_LINKAGE.md` does not already carry it.

PART C - the compile-scan fix
- `00_install.R` flags `lattice`, `Matrix` and `survival` as missing because it searches only the
  project library. They are base-recommended packages installed with R itself.
- Fix it to search the full library path (`.libPaths()`, including the system library) when deciding
  whether a package is present. Work out the correct set of base/recommended packages
  programmatically rather than hard-coding three names — others will hit this.
- Verify: the scan reports nothing missing on this machine, and still correctly reports a genuinely
  absent package. Test both directions.

PART D - the Bioconductor repository URL
- In `renv.lock`, change the Bioconductor repository entry from the `r-universe.dev` URL to
  `https://bioconductor.org/packages/3.22/bioc`, matching the pinned release.
- Confirm all three Bioconductor packages (`BiocVersion`, `edgeR`, `limma`) resolve at EXACTLY their
  pinned versions from the new URL, and that no other package's source changes.
- Show the diff plainly in the PR body — this changes where packages come from, and it was approved
  on the understanding the diff would be reviewed before merging.
- Do NOT change any package VERSION. If a pinned version is unavailable from the new repository,
  STOP and report rather than adjusting the pin.

PART E - PR
- Push `gyd/docs` and open a PR against `origin/main`. Title: "Docs: C1/C2 reports and the prompt
  series; compile-scan and Bioconductor fixes". Body: what was stranded and why; what was brought
  across and what deliberately was not; the dated notes added; the two fixes with the Bioconductor
  diff shown in full; and an explicit statement that no analysis, constant or fitted number changes.
- Do not push to `main`. Verify `origin/main` unchanged.
- Note in the body that `gyd/c1-c2-review` (PR #1) can be closed and the branch deleted once this
  merges — but confirm that by diffing against the merged result first, and leave the deletion to the
  repository owner.

VERIFY (report all)
1. `git diff --stat origin/main...gyd/docs` — confirm only the expected paths.
2. What the C1/C2 reports read from `runs/`, established by parsing; what was tracked and what
   excluded; whether they render from a clean clone.
3. The stale-documentation corrections made, listed individually.
4. The dated notes added, and what each records.
5. Compile scan: nothing falsely missing on this machine; a genuinely absent package still detected.
6. Bioconductor: all three packages resolve at exactly their pinned versions from the new URL; no
   other source changed; the diff.
7. Whether `gyd/c1-c2-review` is fully superseded by this branch.
8. PR opened; `origin/main` unchanged; `results/`, `data/` and `scripts/` (beyond `00_install.R`)
   untouched by checksum.

CONSTRAINTS
- No analysis, no constant, no fitted number, no manuscript edit.
- `runs/` trees do not come across. Only the text tables the reports actually read, if any.
- Correct stale documentation; ANNOTATE superseded findings. Do not rewrite a report to match what we
  learned later — the audit trail is the point.
- Do not change any package version in `renv.lock`.
- Autonomous; commit in parts: "docs: bring the C1/C2 reports and prompt series onto main",
  "docs: refresh indices and runbook for the current layout",
  "docs: note where C1/C2 findings have since been superseded",
  "install: search the full library path in the compile scan",
  "renv: point the Bioconductor repository at bioconductor.org".
```
