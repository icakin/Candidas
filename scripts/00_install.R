#!/usr/bin/env Rscript
# =============================================================================
# 00_install.R - install and verify the R side of the Candidas environment
# =============================================================================
# IDEMPOTENT. Safe to run repeatedly; it only installs what is missing.
#
#   Rscript scripts/00_install.R            # install + verify (compiles a toy Stan model)
#   Rscript scripts/00_install.R --no-stan  # skip the Stan compile check (fast)
#
# What it does, in order:
#   0. CHECKS SYSTEM PREREQUISITES before touching the library, and stops with
#      the specific missing item: Xcode tools -> Fortran toolchain -> Homebrew
#      system libraries -> conda interference -> Bioconductor configuration.
#      Checking first matters because the alternative is a build failing deep in
#      a dependency tree with a linker error nobody can read.
#   1. checks the R version against renv.lock
#   2. installs renv itself if absent
#   3. configures Bioconductor (BiocManager, the pinned release) BEFORE the
#      restore, so edgeR/limma can be found at all
#   4. restores the project library from renv.lock (this is the pin - it covers
#      CRAN *and* Bioconductor packages, edgeR included)
#   5. verifies every package the eighteen numbered scripts actually load
#   6. detects which Stan backend brms will use, then COMPILES AND SAMPLES a
#      10-line toy model to prove the C++ toolchain works. Fails loudly with
#      instructions if it cannot.
#
#   Rscript scripts/00_install.R --no-stan     skip the Stan compile check
#   Rscript scripts/00_install.R --check-only  run step 0 and stop
#   Rscript scripts/00_install.R --force       report prerequisite problems but
#                                              attempt the restore anyway
#
# It does NOT touch the Python side (the etc-GEM engine in cauris_etcgem/ and
# matplotlib for 17_schematic.py). See SETUP.md.
# =============================================================================

args       <- commandArgs(trailingOnly = TRUE)
SKIP_STAN  <- "--no-stan"    %in% args
CHECK_ONLY <- "--check-only" %in% args
FORCE      <- "--force"      %in% args

# ---- how packages are OBTAINED (not which versions) -------------------------
# renv.lock is the pin and remains the sole authority on WHICH version of each
# package is installed. Nothing below can change that: these settings only
# affect where a given version is fetched from and whether a prebuilt binary is
# preferred over the source tarball of the SAME version.
#
# Preferring binaries is worth doing because most of what needs a compiler here
# needs it only because CRAN serves binaries for the CURRENT version alone, so
# every superseded pin falls back to source. It is a partial remedy, not a cure
# (see the note printed by check_binary_coverage() below).
if (Sys.info()[["sysname"]] == "Darwin") {
  # respect a deliberate user setting; only fill in the default
  if (is.null(getOption("pkgType")) || identical(getOption("pkgType"), "source"))
    options(pkgType = "binary")
  if (is.null(getOption("install.packages.compile.from.source")))
    options(install.packages.compile.from.source = "never")
}
options(repos = c(CRAN = "https://cloud.r-project.org"))

.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  d  <- if (length(fa)) dirname(normalizePath(
    # R replaces every space in --file= with "~+~", so a project path that
    # contains a space comes back mangled. Un-mangle it before normalising.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE),
    mustWork = FALSE))
        else tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (length(d) == 0 || is.na(d) || !nzchar(d)) d <- getwd()
  normalizePath(d, mustWork = FALSE)
})
PROJ <- dirname(.this_dir)
LOCK <- file.path(PROJ, "renv.lock")

hdr <- function(x) message("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74))
ok  <- function(...) message("  [ok]   ", ...)
bad <- function(...) message("  [FAIL] ", ...)

hdr("00_install.R - Candidas R environment")
message("  project : ", PROJ)
message("  R       : ", getRversion(), "  (", R.version$platform, ")")


# =============================================================================
# 0. SYSTEM PREREQUISITES  (macOS; see SETUP.md "System prerequisites")
# =============================================================================
# Ordered so the cheapest and most fundamental check comes first, and so the
# first thing that fails is the first thing to fix. Each failure names the
# missing item and how to install it.
#
# NOT VERIFIED ON A CLEAN MACHINE. These checks encode a reported clean-Mac
# failure plus a dependency analysis of renv.lock; they have not been confirmed
# against a machine that actually lacks the items. Use --force to proceed anyway.

.PREREQ_PROBLEMS <- character(0)
prob <- function(what, why, fix) {
  bad(what)
  message("           ", why)
  message("           fix: ", fix)
  .PREREQ_PROBLEMS <<- c(.PREREQ_PROBLEMS, what)
}
.has_cmd <- function(x) nzchar(Sys.which(x))

if (identical(Sys.info()[["sysname"]], "Darwin")) {
  hdr("0/6  system prerequisites (macOS)")

  # ---- 0.1 Xcode command line tools ----------------------------------------
  clt <- suppressWarnings(system2("xcode-select", "-p", stdout = TRUE, stderr = FALSE))
  if (length(clt) && nzchar(clt[1]) && dir.exists(clt[1])) {
    ok("Xcode command line tools  (", clt[1], ")")
  } else {
    prob("Xcode command line tools missing",
         "clang / clang++ / make are required by every source build, and by Stan.",
         "xcode-select --install")
  }

  # ---- 0.2 the CRAN gfortran toolchain -------------------------------------
  # Six pinned packages link the Fortran runtime and have no binary at their pin
  # (Matrix, mgcv, nlme, mvtnorm, nleqslv, edgeR), so they build from source.
  fc <- tryCatch(system("R CMD config FC", intern = TRUE, ignore.stderr = TRUE)[1],
                 error = function(e) "")
  fc_bin <- if (length(fc) && nzchar(fc)) sub(" .*$", "", fc) else ""
  fc_ok  <- nzchar(fc_bin) && (file.exists(fc_bin) || .has_cmd(basename(fc_bin)))
  if (fc_ok) {
    ok("Fortran compiler         (", fc_bin, ")")
    if (grepl("^/opt/homebrew|^/usr/local/Cellar", fc_bin))
      message("  [warn] that looks like Homebrew's gfortran. CRAN-built R expects the",
              "\n         CRAN toolchain from https://mac.r-project.org/tools/ ; a Homebrew",
              "\n         gfortran can link but produce packages that misbehave.")
  } else {
    prob("Fortran compiler not found",
         paste0("`R CMD config FC` gives '", fc, "'. Matrix, mgcv, nlme, mvtnorm, ",
                "nleqslv and edgeR\n           build from source at their pinned versions ",
                "and all link libgfortran."),
         "install the CRAN GNU Fortran package from https://mac.r-project.org/tools/\n                (NOT Homebrew's gcc - see SETUP.md)")
  }

  # ---- 0.3 Homebrew system libraries ---------------------------------------
  # Only checked, never installed: installing system packages is the user's call.
  #
  # Names are matched allowing a version suffix, because Homebrew ships several
  # of these as versioned, keg-only formulae - openssl is `openssl@3`, icu4c is
  # `icu4c@78`. Matching the bare name exactly reports them missing when they
  # are present, which is worse than not checking at all.
  #
  # REQUIRED are the two that were actually reported as failing. The rest are
  # ADVISORY: macOS ships libxml2, and stringi will bundle its own ICU (slowly)
  # rather than fail, so a missing one of those is a slow build, not a dead one.
  brew_required <- c(
    openssl  = "openssl (2.3.5), curl (7.0.0)",
    freetype = "ragg (1.5.0), textshaping (1.0.4), systemfonts (1.3.1)")
  brew_advised <- c(
    harfbuzz     = "textshaping",
    fribidi      = "textshaping",
    fontconfig   = "systemfonts",
    libpng       = "ragg",
    "jpeg-turbo" = "ragg",
    libtiff      = "ragg",
    libxml2      = "xml2 (macOS also ships one)",
    icu4c        = "stringi (bundles its own if absent)")

  if (.has_cmd("brew")) {
    installed <- tryCatch(system2("brew", c("list", "--formula", "-1"),
                                  stdout = TRUE, stderr = FALSE),
                          error = function(e) character(0))
    # present if any installed formula is `name` or `name@<version>`
    have_brew <- function(n) any(grepl(paste0("^", n, "(@.*)?$"), installed))

    miss_req <- names(brew_required)[!vapply(names(brew_required), have_brew, logical(1))]
    miss_adv <- names(brew_advised)[!vapply(names(brew_advised), have_brew, logical(1))]

    if (!length(miss_req)) {
      ok("Homebrew: required libraries present  (openssl, freetype)")
    } else {
      prob(paste0("Homebrew libraries missing: ", paste(miss_req, collapse = ", ")),
           paste0("needed by: ",
                  paste(sprintf("%s -> %s", miss_req, brew_required[miss_req]), collapse = "; ")),
           paste0("brew install ", paste(miss_req, collapse = " ")))
    }
    if (length(miss_adv))
      message("  [note] optional Homebrew libraries absent: ",
              paste(miss_adv, collapse = ", "),
              "\n         (", paste(sprintf("%s -> %s", miss_adv, brew_advised[miss_adv]),
                                    collapse = "; "), ")",
              "\n         These make a source build slower or fall back to a bundled copy;",
              "\n         they do not stop the restore. brew install ",
              paste(miss_adv, collapse = " "))
  } else {
    prob("Homebrew not found",
         paste0("openssl and freetype are needed by packages that build from source at\n",
                "           their pinned versions: openssl, curl, ragg, textshaping, systemfonts."),
         "install Homebrew from https://brew.sh, then:\n                brew install openssl freetype")
  }

  # ---- 0.4 conda on PATH ahead of the system libraries ---------------------
  # The reported curl/libkrb5 failure. Detect it rather than let the build die.
  conda_pref <- nzchar(Sys.getenv("CONDA_PREFIX"))
  path_parts <- strsplit(Sys.getenv("PATH"), ":", fixed = TRUE)[[1]]
  conda_on_path <- grep("(conda|miniforge|mambaforge|anaconda)", path_parts,
                        value = TRUE, ignore.case = TRUE)
  if (conda_pref || length(conda_on_path)) {
    prob("conda is active on PATH",
         paste0("conda's lib directory ahead of the system libraries makes `curl` compile\n",
                "           against conda's libkrb5 and link against the system libcurl. ",
                "That is\n           the reported clean-Mac failure. CONDA_PREFIX='",
                Sys.getenv("CONDA_PREFIX"), "'"),
         "conda deactivate   (repeat until no (env) prefix), then re-run;\n                or: env PATH=\"/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin\" Rscript scripts/00_install.R")
  } else {
    ok("conda not on PATH        (no libkrb5 interference expected)")
  }
} else {
  message("\n  Not macOS (", Sys.info()[["sysname"]], ") - the prerequisite checks are ",
          "macOS-specific and are skipped.\n  Linux: build-essential. Windows: Rtools. ",
          "Neither is documented or tested for this project.")
}

# ---- 0.5 how much of the lockfile can come from a binary --------------------
# Informational, and the honest counterweight to options(pkgType = "binary"):
# a binary preference cannot help where no binary of the PINNED version exists.
check_binary_coverage <- function() {
  ap <- tryCatch(available.packages(type = "binary"), error = function(e) NULL)
  if (is.null(ap) || !nrow(ap)) { message("  [warn] could not reach the CRAN binary index"); return(invisible()) }
  lk <- tryCatch(jsonlite::fromJSON(LOCK)$Packages, error = function(e) NULL)
  if (is.null(lk)) return(invisible())
  want  <- vapply(lk, function(p) p$Version, character(1))
  exact <- vapply(names(want), function(n)
    n %in% rownames(ap) && identical(unname(ap[n, "Version"]), unname(want[[n]])), logical(1))
  message(sprintf("  binary coverage of the pin: %d of %d packages have an exact CRAN binary;",
                  sum(exact), length(want)))
  message(sprintf("  %d must build from source. That is why the toolchain above is required.",
                  sum(!exact)))
}
if (requireNamespace("jsonlite", quietly = TRUE)) check_binary_coverage()

if (length(.PREREQ_PROBLEMS)) {
  message("\n  ", strrep("-", 70))
  message("  ", length(.PREREQ_PROBLEMS), " prerequisite problem(s) found:")
  for (p in .PREREQ_PROBLEMS) message("    - ", p)
  message("  ", strrep("-", 70))
  message("  Full explanation: SETUP.md, 'System prerequisites'.")
  if (!FORCE)
    stop("Prerequisites not met. Fix the above, or re-run with --force to try anyway.",
         call. = FALSE)
  message("  --force given: continuing despite the above.")
} else if (identical(Sys.info()[["sysname"]], "Darwin")) {
  ok("all system prerequisites present")
}

if (CHECK_ONLY) {
  hdr("--check-only: stopping before the restore")
  quit(status = if (length(.PREREQ_PROBLEMS)) 1L else 0L)
}


# ---- 1. R version vs the lockfile -------------------------------------------

if (file.exists(LOCK)) {
  lock_r <- tryCatch({
    j <- readLines(LOCK, warn = FALSE)
    i <- grep('"Version"', j)[1]
    sub('.*"Version"\\s*:\\s*"([^"]+)".*', "\\1", j[i])
  }, error = function(e) NA_character_)
  if (!is.na(lock_r)) {
    if (identical(as.character(getRversion()), lock_r)) {
      ok("R version matches renv.lock (", lock_r, ")")
    } else {
      message("  [warn] renv.lock was built on R ", lock_r, "; you are on ",
              getRversion(), ".\n",
              "         renv will still restore, but binaries may be rebuilt from source\n",
              "         and exact numeric reproduction is no longer guaranteed.")
    }
  }
} else {
  bad("renv.lock not found at ", LOCK)
  stop("Cannot pin the environment without renv.lock. Is this the project root?")
}


# ---- 2. renv ----------------------------------------------------------------

hdr("1/6  renv")
if (!requireNamespace("renv", quietly = TRUE)) {
  message("  installing renv ...")
  install.packages("renv")
}
if (!requireNamespace("renv", quietly = TRUE))
  stop("renv could not be installed. Check your internet connection and CRAN mirror.")
ok("renv ", as.character(packageVersion("renv")))


# ---- 2b. Bioconductor, BEFORE the restore -----------------------------------
# edgeR (needed by 13_capacity_expression.R) and limma are Bioconductor, not
# CRAN. Without the Bioconductor repositories configured, renv cannot FIND them
# at all - a different failure from "the compiler is missing", and one no
# toolchain fixes. BiocManager must therefore exist and the release be set
# BEFORE restore() runs.

hdr("2/6  Bioconductor  (edgeR, limma - needed by 13_capacity_expression.R)")

BIOC_VERSION <- tryCatch({
  j <- jsonlite::fromJSON(LOCK)
  if (!is.null(j$Bioconductor$Version)) j$Bioconductor$Version else NA_character_
}, error = function(e) NA_character_)

if (is.na(BIOC_VERSION)) {
  message("  [warn] renv.lock records no Bioconductor version. edgeR may resolve to",
          "\n         whatever release BiocManager defaults to, which is not a pin.",
          "\n         Report this rather than working around it.")
} else {
  ok("Bioconductor release pinned in renv.lock: ", BIOC_VERSION)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("  installing BiocManager ...")
  install.packages("BiocManager")
}
if (requireNamespace("BiocManager", quietly = TRUE)) {
  ok("BiocManager ", as.character(packageVersion("BiocManager")))
  if (!is.na(BIOC_VERSION)) {
    # Prefer the CANONICAL Bioconductor mirror for the pinned release.
    #
    # WHY: renv.lock records these packages' Repository as
    # https://bioc-release.r-universe.dev, which serves none of BiocVersion
    # 3.22.0 / edgeR 4.8.2 / limma 3.66.0 as a macOS arm64 binary, so they build
    # from source (which is why edgeR needs gfortran). The canonical mirror
    # https://bioconductor.org/packages/<release>/bioc has all three at exactly
    # the pinned versions as arm64 binaries.
    #
    # This changes WHERE a package comes from, never WHICH version: renv.lock
    # still decides that, and a mirror that lacks the pinned version is simply
    # skipped. Correcting the Repository field in renv.lock would be the proper
    # fix and belongs in its own change.
    bioc_repos <- c(
      BioCsoft  = sprintf("https://bioconductor.org/packages/%s/bioc", BIOC_VERSION),
      BioCann   = sprintf("https://bioconductor.org/packages/%s/data/annotation", BIOC_VERSION),
      BioCexp   = sprintf("https://bioconductor.org/packages/%s/data/experiment", BIOC_VERSION))
    options(repos = c(getOption("repos"), bioc_repos))
    ok("Bioconductor repositories configured for release ", BIOC_VERSION)
    have <- tryCatch(available.packages(repos = bioc_repos[["BioCsoft"]], type = "binary"),
                     error = function(e) NULL)
    if (!is.null(have) && "edgeR" %in% rownames(have))
      message("    edgeR available there: ", have["edgeR", "Version"],
              "  (binary, ", R.version$platform, ")")
  }
} else {
  bad("BiocManager could not be installed - edgeR and limma will not restore.")
}


# ---- 3. restore the project library -----------------------------------------

hdr("3/6  renv::restore()  (installs every pinned package, CRAN + Bioconductor)")
message("  On a fresh machine this is the slow step. rstan and StanHeaders DO have")
message("  binaries at their pinned versions, so the worst of it is avoided; what")
message("  still compiles is the ~73 pins CRAN no longer serves as binaries")
message("  (RcppParallel builds a bundled TBB and is the slowest of them).")
message("  It is idempotent - already-correct packages are skipped.\n")

restore_ok <- tryCatch({
  renv::restore(project = PROJ, prompt = FALSE)
  TRUE
}, error = function(e) { bad("renv::restore() failed: ", conditionMessage(e)); FALSE })

if (!restore_ok) {
  message("\n  If restore failed on a Bioconductor package (edgeR / limma), install")
  message("  BiocManager first and retry:")
  message("      install.packages('BiocManager'); BiocManager::install(version = '3.22')")
}


# ---- 4. can every package the pipeline loads actually be loaded? ------------

hdr("4/6  package check  (everything scripts/*.R actually loads)")

# Kept in sync by hand with
#   grep -hoE '(library|require|requireNamespace)\(|\b[A-Za-z][A-Za-z0-9._]*::' scripts/*.R
NEEDED <- c(
  # tidyverse core + plotting
  "tidyverse", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr",
  "ggplot2", "scales", "ggrepel", "ggdist", "patchwork", "rlang", "grid",
  # fitting / numerics
  "minpack.lm", "zoo",
  # the four click-driven apps (01, 04, 05, 06)
  "shiny", "rstudioapi",
  # Bayesian stage (09, 10, 11, 18)
  "brms", "posterior",
  # 13_capacity_expression.R (Bioconductor)
  "edgeR",
  # xlsx ingest in 01, and the measured T_internal trace read by 08
  "readxl"
)

status <- vapply(NEEDED, function(p)
  isTRUE(suppressWarnings(requireNamespace(p, quietly = TRUE))), logical(1))

for (p in NEEDED) {
  if (status[[p]]) ok(sprintf("%-12s %s", p, as.character(packageVersion(p))))
  else             bad(sprintf("%-12s MISSING", p))
}

if (any(!status)) {
  missing <- names(status)[!status]
  message("")
  bad(length(missing), " package(s) missing: ", paste(missing, collapse = ", "))
  if ("edgeR" %in% missing)
    message("      edgeR is Bioconductor, not CRAN:\n",
            "          install.packages('BiocManager')\n",
            "          BiocManager::install('edgeR', version = '3.22')")
  stop("Environment incomplete - fix the above before running the pipeline.")
}
ok("all ", length(NEEDED), " packages present")


# ---- 5. Stan backend + toolchain --------------------------------------------

hdr("5/6  Stan backend")

# WHICH BACKEND? 09_bayesian_models.R calls brms::brm() without a `backend`
# argument and nothing in the project sets options(brms.backend), so brms uses
# its own default. Report what that resolves to on THIS machine.
opt_backend <- getOption("brms.backend", NULL)
has_rstan   <- requireNamespace("rstan",    quietly = TRUE)
has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE)

backend <- if (!is.null(opt_backend)) as.character(opt_backend) else "rstan"

message("  options(brms.backend) : ", if (is.null(opt_backend)) "<unset> -> brms default" else opt_backend)
message("  rstan installed       : ", has_rstan,
        if (has_rstan) paste0("  (", packageVersion("rstan"), ")") else "")
message("  cmdstanr installed    : ", has_cmdstan,
        if (has_cmdstan) paste0("  (", packageVersion("cmdstanr"), ")") else "")
message("  ==> 09_bayesian_models.R will use backend: ", backend)

if (identical(backend, "rstan") && !has_rstan)
  stop("brms will use the rstan backend but rstan is not installed.\n",
       "    install.packages('rstan', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))")

if (has_rstan) {
  message("  StanHeaders           : ", as.character(packageVersion("StanHeaders")))
  message("  Stan library          : ", tryCatch(rstan::stan_version(), error = function(e) "?"))
}
message("  C++ compiler          : ", tryCatch(
  system("R CMD config CXX17", intern = TRUE)[1], error = function(e) "?"))
message("  CXX17FLAGS            : ", tryCatch(
  system("R CMD config CXX17FLAGS", intern = TRUE)[1], error = function(e) "?"))

if (SKIP_STAN) {
  message("\n  --no-stan given: skipping the compile-and-sample check.")
} else {
  message("\n  Compiling and sampling a 10-line toy model (first run: 30-90 s) ...")

  toy <- "
data { int<lower=0> N; vector[N] y; }
parameters { real mu; real<lower=0> sigma; }
model { y ~ normal(mu, sigma); }
"
  set.seed(1)
  ydat <- rnorm(50, 3, 1)

  res <- tryCatch({
    if (identical(backend, "cmdstanr") && has_cmdstan) {
      f  <- tempfile(fileext = ".stan"); writeLines(toy, f)
      md <- cmdstanr::cmdstan_model(f)
      fit <- md$sample(data = list(N = 50L, y = ydat), chains = 1,
                       iter_warmup = 100, iter_sampling = 100,
                       seed = 1, refresh = 0, show_messages = FALSE)
      s <- fit$summary(c("mu", "sigma"))
      c(mu = s$mean[1], sigma = s$mean[2])
    } else {
      fit <- rstan::stan(model_code = toy, data = list(N = 50L, y = ydat),
                         chains = 1, iter = 200, refresh = 0, seed = 1)
      m <- rstan::summary(fit)$summary[, "mean"]
      c(mu = unname(m["mu"]), sigma = unname(m["sigma"]))
    }
  }, error = function(e) e)

  if (inherits(res, "error")) {
    bad("Stan could not compile or sample.")
    message("\n  ERROR: ", conditionMessage(res))
    message("\n  ", strrep("-", 70))
    message("  THE BAYESIAN STAGE (09_bayesian_models.R) CANNOT RUN. Fix this first.")
    message("  ", strrep("-", 70))
    message("  macOS   : install the Xcode command-line tools ->  xcode-select --install")
    message("            then, if rstan still fails, follow")
    message("            https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac")
    message("  Linux   : sudo apt install build-essential   (needs g++ and GNU make)")
    message("  Windows : install Rtools matching your R version")
    message("  Or switch to cmdstanr:")
    message("      install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))")
    message("      cmdstanr::install_cmdstan()")
    message("      options(brms.backend = 'cmdstanr')   # in ~/.Rprofile or before sourcing 09")
    stop("Stan backend verification FAILED.")
  }

  ok(sprintf("toy model sampled: mu = %.4f (truth 3), sigma = %.4f (truth 1)",
             res[["mu"]], res[["sigma"]]))
}


# ---- done -------------------------------------------------------------------

hdr("R environment ready")
message("  Next:")
message("    1. Python side :  python3 -m pip install numpy matplotlib   (for 17_schematic.py)")
message("    2. Quarto PDF  :  quarto install tinytex")
message("    3. Run it all  :  bash scripts/run_all.sh")
message("")
message("  See SETUP.md for prerequisites, expected runtimes and disk.")
invisible(TRUE)
