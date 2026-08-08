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
#   1. checks the R version against renv.lock
#   2. installs renv itself if absent
#   3. restores the project library from renv.lock (this is the pin - it covers
#      CRAN *and* Bioconductor packages, edgeR included)
#   4. verifies every package the eighteen numbered scripts actually load
#   5. detects which Stan backend brms will use, then COMPILES AND SAMPLES a
#      10-line toy model to prove the C++ toolchain works. Fails loudly with
#      instructions if it cannot.
#
# It does NOT touch the Python side (the etc-GEM engine in cauris_etcgem/ and
# matplotlib for 17_schematic.py). See SETUP.md.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

args       <- commandArgs(trailingOnly = TRUE)
SKIP_STAN  <- "--no-stan" %in% args

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

hdr("1/4  renv")
if (!requireNamespace("renv", quietly = TRUE)) {
  message("  installing renv ...")
  install.packages("renv")
}
if (!requireNamespace("renv", quietly = TRUE))
  stop("renv could not be installed. Check your internet connection and CRAN mirror.")
ok("renv ", as.character(packageVersion("renv")))


# ---- 3. restore the project library -----------------------------------------

hdr("2/4  renv::restore()  (installs every pinned package, CRAN + Bioconductor)")
message("  On a fresh machine this is the slow step (10-40 min: rstan and StanHeaders")
message("  compile). It is idempotent - already-correct packages are skipped.\n")

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

hdr("3/4  package check  (everything scripts/*.R actually loads)")

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

hdr("4/4  Stan backend")

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
