# Claude Code prompt — Handover: write the handover note and push Candidas (and its submodule) for review (autonomous)

Run from the project root (`.../Candidas TPC/Candidas`). Handover and git only. It runs no analysis,
refits nothing, renders nothing new and changes no result. The purpose is to leave the repository in a
state its author can pick up cold, and to get the C1/C2/C2b work onto GitHub for review.

NOTE TO USER: launch in an auto-approving mode. Fast — a written note and git. The push steps need
network and credentials for `git@github.com:icakin/Candidas.git` AND
`git@github.com:icakin/cauris_etcgem.git`.

CONTEXT. C1, C2 and C2b are complete; local `main` is 22 commits ahead of `origin/main`. C2b already
tidied the tree (run trees consolidated under `runs/` with a manifest, `prompts/` tracked, reports in
`reports/`), so this prompt does NOT repeat that work.

THIS IS SOMEONE ELSE'S REPOSITORY, and so is the submodule. Push BRANCHES and open pull requests. Do
not push to `origin/main` or to the submodule's `origin/master`.

TWO GIT COMPLICATIONS, both of which will bite if handled in the wrong order:

1. **The submodule has its own local commits.** `cauris_etcgem` is at `c224c3c` on branch `master`
   with at least three unpushed commits from C1 and C2b (`0134dd1` gitignore the DE fit cache,
   `04c6107` move the C1 etc-GEM outputs under `runs/C1/`, `c224c3c` remove a tracked .DS_Store).
   The submodule MUST be pushed FIRST. If the superproject is pushed first it records a submodule
   commit that does not exist on the remote, and anyone cloning gets a dangling pointer.
2. **`.gitmodules` declares an https URL while the working remote is ssh**
   (`https://github.com/icakin/cauris_etcgem.git` vs `git@github.com:icakin/cauris_etcgem.git`), and
   the repository is PRIVATE — an unauthenticated `git clone --recursive` fails on it. That is the
   reproducibility blocker C1 documented, and it needs stating in the handover note rather than
   silently fixing.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: prompts/README.md,
reports/{REPRODUCTION.md,N0_SENSITIVITY.md,RUNBOOK.md}, reports/n0/ (the rendered N0 report),
runs/MANIFEST.md, SETUP.md, .gitmodules, and both `git log --oneline origin/main..HEAD` and the
equivalent inside `cauris_etcgem`. Do NOT re-run any analysis, re-render any report, or alter any
number.

PART A - pre-flight
- Confirm the working tree is clean in BOTH the superproject and the submodule, and that C2b's tidy
  landed (runs/ consolidated with a manifest, prompts/ tracked, no stray Rplots.pdf or .DS_Store).
  If anything is uncommitted, list it and stop rather than committing it blind.
- Confirm `results/`, `data/`, `manuscript/`, `outputs/supp_data/` and
  `outputs/supp_data_original_backup/` are byte-identical to their state at the last C2b commit —
  checksum, not assertion.
- Report the exact commit ranges to be pushed, for both repositories.

PART B - push the submodule FIRST
- In `cauris_etcgem`: create a branch off its current `master` named `gyd/c1-c2-review` and push it to
  `origin`. Do NOT push to `origin/master`.
- Open a PR against `origin/master` if `gh` is available; otherwise print the URL, title and body to
  paste. The body should say what changed (C1 output trees moved under `runs/C1/`, gitignore for the
  DE fit cache, .DS_Store removal) and, importantly, that the etc-GEM outputs in `outputs/supp_data/`
  were found NOT to reproduce — capacity regenerates 1.33-1.37x higher, traced to a change in the
  `ANCHOR_MU` bisection, and `figure3_data/model_curves.csv` uses a temperature grid the current code
  cannot emit. Point at `reports/REPRODUCTION.md` in the superproject for the detail.
- Verify the branch landed before going on. If the push fails, STOP — do not push the superproject.

PART C - then push the superproject
- Create a branch off current `main` named `gyd/c1-c2-review` and push it to `origin`.
- Verify the submodule commit the branch records now EXISTS on the submodule's remote. State how you
  verified it.
- Open a PR against `origin/main` (or print URL/title/body). The body should be a condensed version of
  HANDOVER.md's status and open-questions sections, and must say plainly that it is for review rather
  than to be merged unexamined.
- Do NOT force-push anything. `origin/main` and the submodule's `origin/master` must be untouched.

PART D - the handover note
Write `HANDOVER.md` at the project root, for a reader returning cold. Under two pages; link out rather
than restate. It must contain:
- **Status**: C1, C2 and C2b each in two or three sentences, with paths to the full detail
  (`reports/REPRODUCTION.md`, `reports/n0/`, `runs/MANIFEST.md`).
- **What is solid, stated FIRST**: the pipeline runs unattended in 49 min from committed raw data; the
  respirometry half reproduces (worst disagreement 2.4e-10, every Bayesian number within 0.4%
  relative); growth-side quantities are invariant across all N0 treatments at exactly 0.000e+00; and
  the headline result — CUE optimum below body temperature in every taxon — holds and if anything
  strengthens.
- **Manuscript status**, as a table, claim by claim: SURVIVES / WEAKENS / FAILS with the number
  behind each. At minimum: growth turns over while respiration does not; CUE optimum below 37 C;
  growth ~2x more temperature-sensitive than respiration; all four C. auris clades cheaper at fever
  than C. parapsilosis; "Clade IV credibly the cheapest" (FAILS — Clade III vs IV becomes 1.00
  [0.920, 1.319]); rho = -0.90 (circular); the capacity correlations (Fig 4e p = 0.54).
- **What is open**, each with its next action: the N0 back-calculation (see the OxygenModel work —
  this is a method-level question, not a Candidas one); `T_internal`, recorded by the PreSens exports
  and stripped by `01_convert_xlsx.R`; the 88 manual exclusions and para having no data above 40 C;
  the etc-GEM half, where nothing absolute reproduces and the decision on whether it stays in the
  paper is deferred until T_internal and C3.
- **The submodule problem**, explicitly: `cauris_etcgem` is private and `.gitmodules` points at an
  https URL that unauthenticated clones cannot use, so `git clone --recursive` fails for anyone
  outside the group. Note what would fix it (make it public, or vendor the outputs) without doing it.
- **Decisions already taken**, so they are not relitigated: the back-calculation is a design choice
  forced by the sealed vial; ARM 3 is a lower bound, not a candidate; the etc-GEM decision is
  deferred; the estimator work moves to the OxygenModel repo.
- **How to resume**: the exact commands, what to read first, and that `prompts/README.md` lists C3-C8
  as described but NOT YET WRITTEN.
- Cross-reference the OxygenModel repo and its `HANDOVER.md`, since the two projects now share the
  same open question.

VERIFY (report all)
1. Both working trees clean; the commit ranges pushed, per repository.
2. Submodule branch pushed BEFORE the superproject; how you verified the recorded submodule commit
   exists on its remote.
3. Both PRs opened, or both URLs/titles/bodies printed.
4. `origin/main` and the submodule's `origin/master` unchanged, verified.
5. `results/`, `data/`, `manuscript/` and the two etc-GEM output directories byte-identical.
6. `HANDOVER.md` exists and covers all seven required sections; state its length.
7. `prompts/README.md` statuses correct (C1, C2, C2b marked run; C3-C8 pending and unwritten).

CONSTRAINTS
- No analysis, no refits, no re-rendering, no number changes.
- Branches and PRs only. Nothing force-pushed. Neither default branch touched.
- Submodule first, superproject second. If the submodule push fails, stop.
- Do not "fix" the `.gitmodules` URL or the private-repo problem — document them. Changing the URL
  without knowing whether the repo will be made public would make the handover note wrong.
- HANDOVER.md must lead with what is solid before what is open. The reader built this.
- Autonomous; commit in parts: "handover: HANDOVER.md for the Candidas C-series",
  "handover: push cauris_etcgem gyd/c1-c2-review", "handover: push gyd/c1-c2-review for review".
```
