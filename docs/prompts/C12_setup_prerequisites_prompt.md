# Claude Code prompt — C12: make renv::restore() work on a clean Mac (Candidas) (autonomous)

Run from the project root (`.../Candidas TPC/Candidas`). Documentation and installer only. It runs no
analysis, changes no constant and alters no number. It fixes a real gap in the packaging we shipped:
`renv::restore()` does not currently work on a clean machine, so the pinned environment we asked
people to use is not actually usable by them.

NOTE TO USER: launch in an auto-approving mode. Small — `SETUP.md`, `scripts/00_install.R`, and a PR.
The packaging has merged, so branch off `main`. Run the OxygenModel companion prompt FIRST and finish
it — its PART A findings on the shared packages will save work here, though this repository's
dependencies are materially different and PART A must still be redone.

CONTEXT. The packaging PR pinned the R environment with renv and told users to run
`renv::restore()`. On a clean Mac it fails: packages attempt to build from source, requiring the CRAN
gfortran toolchain and Homebrew's `openssl` and `freetype`, and `curl` fails to build when a conda
installation puts its own `libkrb5` ahead of the system libraries on `PATH`. The practical
consequence is that the pipeline has been run against the SYSTEM library rather than the pinned one,
which quietly defeats the purpose of pinning.

TWO THINGS MAKE THIS REPOSITORY HARDER THAN OxygenModel, and they are the reason this is a separate
prompt rather than a re-run:

  * `brms` requires a working Stan backend. That is the heaviest compile in either project, it is the
    most likely single point of failure on a clean machine, and `09_bayesian_models.R` cannot run
    without it. Establish which backend the pinned setup uses (rstan or cmdstanr) — they have quite
    different installation stories, and cmdstanr is not on CRAN at all.
  * `edgeR` comes from BIOCONDUCTOR, not CRAN. That is a different failure mode from compilation:
    `renv::restore()` needs the Bioconductor repositories configured and the right Bioc release
    pinned, or it will not find the package regardless of the toolchain. `13_capacity_expression.R`
    depends on it.

This is our gap, not the reporter's. WE CANNOT TEST THE FIX HERE — there is no clean Mac available —
so write instructions that are correct, honest about what has and has not been verified, and easy for
someone else to confirm. Do not claim the restore works.

---

```
Work AUTONOMOUSLY; commit in parts; print a summary. Read first: `SETUP.md`, `scripts/00_install.R`,
`renv.lock`, `.Rprofile`, `renv/activate.R`, `renv/settings.json` (for repository configuration), and
the library() calls across all eighteen scripts. Branch off `main`:
`git switch -c gyd/setup-prereqs main`. Do NOT run the pipeline.

PART A - work out what actually needs building, for THIS lockfile
- From `renv.lock`, identify every package with compiled code and classify what each needs: a Fortran
  compiler, a named system library, or nothing beyond the Xcode command line tools. Report the list
  and how you determined it. Do not assume the OxygenModel answer transfers.
- STAN: establish which backend is pinned and how it is obtained. Report its build requirements,
  expected compile time, whether a macOS binary exists for the pinned version, and — if cmdstanr —
  that it is not a CRAN package and how `renv` is expected to restore it.
- BIOCONDUCTOR: establish how `edgeR` and its dependencies are recorded in `renv.lock` (source field,
  repository entry), whether the Bioconductor release is pinned, and what a clean `renv::restore()`
  needs in order to find them. If `BiocManager` is required first, say so and put it in the right
  order in the installer.
- Note which packages have CRAN macOS binaries, since the PART C fallback makes those a non-issue.
  Be explicit that Bioconductor binaries are a separate question.

PART B - SETUP.md
Add a PREREQUISITES section BEFORE the install instructions:
- Xcode command line tools (`xcode-select --install`).
- The CRAN gfortran toolchain for macOS, with the download URL, and a note that Homebrew's gcc is NOT
  a substitute for CRAN-built R.
- Homebrew packages: `openssl` and `freetype` at minimum, plus anything PART A adds, each labelled
  with the R package it serves.
- A STAN subsection: what the backend needs, roughly how long it takes, and how to verify it works
  (compile and sample a trivial model) BEFORE discovering the problem partway through
  `09_bayesian_models.R`.
- A BIOCONDUCTOR subsection: what is needed for `edgeR` to restore, and the release pinned.
- A KNOWN ISSUES subsection covering the conda `libkrb5` problem: symptom (curl build failure), cause
  (conda ahead of system libraries on `PATH`), workaround (deactivate conda, or prepend system paths,
  for the duration of the restore).
- State plainly that these prerequisites were derived from a reported clean-Mac failure and have NOT
  been confirmed on a clean machine here.
- Keep Linux and Windows out of scope unless already documented.

PART C - a fallback that avoids the problem
- Prefer binaries on macOS in `scripts/00_install.R`: `options(pkgType = "binary")` and
  `install.packages.compile.from.source = "never"` where appropriate, guarded to macOS and not
  overriding a deliberate user setting.
- Make `00_install.R` CHECK prerequisites BEFORE attempting a restore, in a sensible order — Xcode
  tools, toolchain, system libraries, Bioconductor configuration, Stan backend — failing with a clear
  message naming the missing item and how to install it.
- Verify the Stan backend compiles by building and sampling a trivial model, and report the result.
  A restore that "succeeds" while Stan cannot compile is worse than one that fails cleanly, because
  the failure then surfaces an hour into a pipeline run.
- The fallback must not change WHICH versions are installed. `renv.lock` is the pin; this changes how
  they are obtained. Say so in a comment.

PART D - the PR
- Push and open a PR against `main`. Title: "SETUP: system prerequisites, Bioconductor and Stan
  checks for renv::restore()". Body: what fails and why; the prerequisite list with the package each
  serves; the Stan and Bioconductor findings; the fallback; an explicit statement that no analysis,
  constant or number changes; and a request that the reporter confirm on a clean machine.
- Do not push to `main`. Verify `main` unchanged.

VERIFY (report all)
1. The compiled-code package list for THIS lockfile, how determined, and how it differs from
   OxygenModel's.
2. The Stan backend: which one, its requirements, whether a pinned-version binary exists, and the
   trivial-model compile result.
3. Bioconductor: how `edgeR` is recorded, whether the release is pinned, and what a clean restore
   needs.
4. The prerequisite list, with the R package each item serves.
5. What `00_install.R` now checks, in what order, and the messages it produces.
6. Confirmation the binary fallback cannot change installed VERSIONS, only their source.
7. That `SETUP.md` states plainly what has not been verified here.
8. PR opened; `main` unchanged.

CONSTRAINTS
- No analysis, no constant, no number, no pipeline run. Documentation and installer only.
- Do not claim the restore has been verified.
- Do not pin new package versions or edit `renv.lock`. If the Bioconductor release is NOT pinned,
  report that as a finding rather than fixing it — pinning it changes what gets installed and belongs
  in its own change.
- Redo PART A for this lockfile. Do not carry the OxygenModel answer across.
- Autonomous; commit in parts: "SETUP: system prerequisites for a clean-machine restore",
  "SETUP: Stan backend and Bioconductor requirements",
  "install: ordered prerequisite checks and binary preference on macOS".
```
