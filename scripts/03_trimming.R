# --------------------------------------------------------------------------
# 03_trimming.R
#
# Oxygen-time-series trimming by main spline curve shape.
# Temperature x isolate design: each curve is identified by (T, OTU, Replicate),
# where OTU is the unique Candida isolate code (1..18). The trimming logic is
# purely curve-shape based and does not depend on the biology of the labels.
#
#   - Removes all observations after MAX_TIME_MIN before trimming and plotting.
#   - Spline spar set to TRIM_SPAR.
#   - Start detection finds the visible oxygen tip before sustained decline.
#   - Endpoint detection uses the FIRST meaningful slope inflection after start.
#   - Fallback endpoint is the global steepest drop after start.
#   - Manual overrides can replace automatic start and/or end.
#
# What this script does:
#   1. Reads tables/Oxygen_All_Long.csv
#   2. Removes Time > MAX_TIME_MIN
#   3. Creates a short code for each curve (C001, C002, ...)
#   4. Automatically trims each curve
#   5. Lets you manually override start and/or end by curve code
#   6. Saves trimmed data and a PDF showing diagnostics
#
# Manual override rule:
#   - If a curve code is NOT entered below, automatic trimming is used.
#   - If a curve code IS entered:
#       * keep_start filled  -> manual start
#       * keep_start = NA    -> automatic start
#       * keep_end filled    -> manual end
#       * keep_end = NA      -> automatic end
#   - Manual keep_end is capped at MAX_TIME_MIN.
# --------------------------------------------------------------------------


# =============================================================================
# 0) Source shared config
# =============================================================================

# Packaging: resolve the script directory from --file= FIRST. Under
# `Rscript scripts/<this>.R` neither rstudioapi nor sys.frame(1)$ofile resolves,
# so this fell back to getwd() and then died on "cannot open file .../config.R".
# Sourcing it from run_all.R was unaffected, which is why the bug stayed hidden.
.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    # R replaces every space in --file= with "~+~", so a project path that
    # contains a space comes back mangled and source() then fails on a path
    # that does not exist. Un-mangle it before normalising.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE),
    mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
})

source(file.path(.this_dir, "config.R"))


# =============================================================================
# 1) Load packages
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(zoo)
  library(stringr)
  library(scales)
})


# =============================================================================
# 2) Main knobs
# =============================================================================

TRIM_SPAR <- 0.2

# Hard cutoff: remove everything after this time before trimming, plotting,
# and saving outputs. The new (temperature x otu) recordings run to ~1450 min
# and the active respiration decline routinely extends past 300 min, so the cap
# is set above the longest recording and the actual fit window is chosen by the
# hybrid start/end logic below rather than by this cap.
MAX_TIME_MIN <- 1455

MIN_ROWS_SERIES <- 8
MIN_RANGE_O2    <- 0.05

# Main descending-region detection
NEG_SLOPE_EPS     <- 0.00015
MIN_RUN_POINTS    <- 6
MAX_GAP_POINTS    <- 2
SEARCH_START_FRAC <- 0.08
SEARCH_END_FRAC   <- 0.98

# Start selection: detect first sustained real decline, then move back to tip
START_BACKOFF_POINTS <- 0
PRETRANS_SLOPE_EPS   <- 0.00005
PRETRANS_MIN_RUN     <- 5
SLOPE_TREND_EPS      <- 0.00001
START_SMOOTH_K       <- 5

# Real-decline start controls
START_LOOKBACK_POINTS <- 200
START_NEG_FRAC        <- 0.60
START_DROP_FRAC       <- 0.02
START_NOT_RISING_FRAC <- 0.70
START_WIGGLE_TOL      <- 0.00005
START_NOISE_MULT      <- 2
START_NOISE_FRAC      <- 0.15

# Endpoint detection
END_SMOOTH_K <- 5

# TRUE = end at first meaningful slope inflection / steepest drop.
# This keeps the first active respiration phase and removes later tails.
END_AT_STEEPEST <- TRUE

# TRUE = prefer the first local minimum of the smoothed slope after start.
# FALSE = use the global steepest slope after start.
END_USE_FIRST_INFLECTION <- TRUE

# Half-window size used to identify a local slope minimum.
END_INFLECTION_WIN <- 5

# Minimum negativity required for a slope inflection to count.
# Increase if small wiggles are chosen too early.
# Decrease if real declines are missed.
END_INFLECTION_NEG_EPS <- 0.0005

# Safety controls to avoid ending immediately after the start.
END_MIN_POINTS_AFTER_START   <- 4
END_MIN_DURATION_AFTER_START <- 6

# Older recovery/plateau-tail behaviour.
# These are only used when END_AT_STEEPEST = FALSE.
END_MIN_GAP    <- 2
END_RECOV_RUN   <- 2
END_RECOV_FRAC  <- 0.99
END_MONO_FRAC   <- 0.40
END_SLOPE_EPS   <- 0.00005

CUT_LONG_PLATEAU_TAIL   <- TRUE
PLATEAU_SLOPE_EPS       <- 0.0005
PLATEAU_RUN_POINTS      <- 8
PLATEAU_MIN_AFTER_START <- 10
PLATEAU_KEEP_POINTS     <- 2

# ---------------------------------------------------------------------------
# Standardised hybrid trimming (preferred path)
# ---------------------------------------------------------------------------
# When END_USE_HYBRID = TRUE the script ignores the older inflection / steepest
# logic and uses a deterministic rule that is reproducible across replicates:
#
#   START: argmax of the smoothed spline within
#          [START_SEARCH_MIN_TIME, START_SEARCH_MAX_TIME] minutes.
#
#   END:   min(
#            start_time + END_MAX_WINDOW_MIN,
#            time at which the smoothed spline first drops by
#            END_DROP_FRAC of (O2(start) - min(O2_fit after start)),
#            end of recording
#          )
#
# Series with fewer than MIN_FIT_POINTS kept points or shorter than
# MIN_FIT_DURATION_MIN minutes are skipped. Manual overrides still take
# precedence if filled in.
# ---------------------------------------------------------------------------

END_USE_HYBRID          <- TRUE

# Start = argmax of the smoothed spline (the oxygen tip). In this dataset the
# tip occurs at roughly 70-140 min, so the search window spans 10-200 min. A
# too-narrow window (the old 25-50 min) pins the start on the rising limb,
# before the real tip.
START_SEARCH_MIN_TIME   <- 10
START_SEARCH_MAX_TIME   <- 200

# End = whichever comes first of:
#   * start + END_MAX_WINDOW_MIN,
#   * the time the decline DECELERATES, i.e. the smoothed slope first recovers to
#     END_SLOPE_RECOVER_FRAC of its steepest (most negative) value after the tip
#     - this stops the window at the onset of the flat floor, which is important
#     for the slow 25 degC curves whose floor is noisy/oscillatory,
#   * the time the spline has fallen by END_DROP_FRAC of its post-tip drop (a
#     safety cap if the slope never clearly recovers),
#   * the end of recording.
END_MAX_WINDOW_MIN      <- 900
END_DROP_FRAC           <- 0.85
# Higher -> the window ends EARLIER (the slope only needs to recover a little
# from its steepest value before the curve is treated as flattening). 0.6 keeps
# the active decline but stops before the floor on the slower/hotter curves.
END_SLOPE_RECOVER_FRAC  <- 0.6

MIN_FIT_DURATION_MIN    <- 20
MIN_FIT_POINTS          <- 15

# PDF axis detail
X_MAJOR_BREAKS_N <- 12
Y_MAJOR_BREAKS_N <- 10

# Diagnostic-plot zoom: the recordings run to ~1450 min but the informative part
# is the tip + active decline. Each diagnostic panel is zoomed to the kept window
# padded by PLOT_PAD_FRAC of that window on each side (at least PLOT_PAD_MIN_MIN
# minutes), so the long flat tail past the chosen end is not shown.
PLOT_PAD_FRAC   <- 0.35
PLOT_PAD_MIN_MIN <- 30


# =============================================================================
# 3) Manual override section
# =============================================================================
# 1. Run the script once.
# 2. Open:
#      tables/Oxygen_Curve_Code_Key.csv
#      figures/oxygen_trimming_diagnostics.pdf
# 3. Find the curve code you want to adjust, e.g. C014.
# 4. Enter the code below and fill only the side(s) to change.
# =============================================================================

manual_trim_overrides <- tibble::tibble(
  curve_code = character(),
  keep_start = numeric(),
  keep_end   = numeric()
)

# Example:
# manual_trim_overrides <- tibble::tibble(
#   curve_code = c("C014", "C021"),
#   keep_start = c(NA, 25),
#   keep_end   = c(90, NA)
# )


# =============================================================================
# 4) Helper functions
# =============================================================================

safe_spline_fit <- function(time, oxygen, spar = TRIM_SPAR) {
  fit <- tryCatch(
    smooth.spline(time, oxygen, spar = spar),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(rep(NA_real_, length(time)))
  }
  
  tryCatch(
    predict(fit, x = time)$y,
    error = function(e) rep(NA_real_, length(time))
  )
}


get_slopes <- function(time, y) {
  dt <- diff(time)
  dy <- diff(y)
  
  s <- dy / dt
  s[!is.finite(s)] <- NA_real_
  
  s
}


smooth_slopes <- function(slopes, k = START_SMOOTH_K) {
  out <- slopes
  
  if (length(slopes) >= k) {
    sm <- zoo::rollmedian(slopes, k = k, fill = NA, align = "center")
    sm[is.na(sm)] <- slopes[is.na(sm)]
    out <- sm
  }
  
  out
}


rle_to_runs <- function(flag_vec) {
  r <- rle(flag_vec)
  
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  
  tibble::tibble(
    value = r$values,
    start = starts,
    end   = ends,
    len   = r$lengths
  )
}


merge_descending_runs <- function(run_tbl, max_gap = MAX_GAP_POINTS) {
  if (nrow(run_tbl) == 0) {
    return(run_tbl)
  }
  
  desc_runs <- run_tbl %>%
    dplyr::filter(value)
  
  if (nrow(desc_runs) <= 1) {
    return(
      desc_runs %>%
        dplyr::select(start, end, len)
    )
  }
  
  out <- list()
  
  cur_start <- desc_runs$start[1]
  cur_end   <- desc_runs$end[1]
  
  for (i in 2:nrow(desc_runs)) {
    gap <- desc_runs$start[i] - cur_end - 1
    
    if (gap <= max_gap) {
      cur_end <- desc_runs$end[i]
    } else {
      out[[length(out) + 1]] <- tibble::tibble(
        start = cur_start,
        end = cur_end
      )
      
      cur_start <- desc_runs$start[i]
      cur_end   <- desc_runs$end[i]
    }
  }
  
  out[[length(out) + 1]] <- tibble::tibble(
    start = cur_start,
    end = cur_end
  )
  
  dplyr::bind_rows(out) %>%
    dplyr::mutate(len = end - start + 1L)
}


score_descending_runs <- function(time, o2_fit, slopes, run_tbl) {
  if (nrow(run_tbl) == 0) {
    return(run_tbl)
  }
  
  run_tbl %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      time_start = time[start],
      time_end   = time[min(end + 1L, length(time))],
      duration   = time_end - time_start,
      y_start    = o2_fit[start],
      y_end      = o2_fit[min(end + 1L, length(o2_fit))],
      total_drop = y_start - y_end,
      mean_slope = mean(slopes[start:end], na.rm = TRUE),
      score      = pmax(duration, 0) * pmax(total_drop, 0)
    ) %>%
    dplyr::ungroup()
}


find_main_descending_run <- function(time, o2_fit) {
  slopes <- get_slopes(time, o2_fit)

  if (length(slopes) < 3 || all(!is.finite(slopes))) {
    return(list(run = NULL, slopes = slopes, all_runs = NULL))
  }

  # Start the descending-run search at each curve's own O2 peak (consumption
  # onset) within [START_SEARCH_MIN_TIME, START_SEARCH_MAX_TIME], so the run
  # start tracks the peak per curve instead of a fixed fraction of the record.
  peak_win <- which(
    is.finite(o2_fit) & is.finite(time) &
      time >= START_SEARCH_MIN_TIME & time <= START_SEARCH_MAX_TIME
  )
  idx_peak <- if (length(peak_win) > 0) {
    peak_win[which.max(o2_fit[peak_win])]
  } else {
    max(1L, floor(length(slopes) * SEARCH_START_FRAC))
  }

  lo <- max(1L, min(idx_peak, length(slopes)))
  hi <- min(length(slopes), ceiling(length(slopes) * SEARCH_END_FRAC))

  neg_flag <- rep(FALSE, length(slopes))
  
  idx_search <- lo:hi
  
  neg_flag[idx_search] <- is.finite(slopes[idx_search]) &
    slopes[idx_search] < -NEG_SLOPE_EPS
  
  run_tbl <- rle_to_runs(neg_flag)
  run_tbl <- merge_descending_runs(run_tbl, max_gap = MAX_GAP_POINTS)
  
  if (nrow(run_tbl) == 0) {
    return(list(run = NULL, slopes = slopes, all_runs = NULL))
  }
  
  run_tbl <- run_tbl %>%
    dplyr::filter(len >= MIN_RUN_POINTS)
  
  if (nrow(run_tbl) == 0) {
    return(list(run = NULL, slopes = slopes, all_runs = NULL))
  }
  
  scored <- score_descending_runs(time, o2_fit, slopes, run_tbl) %>%
    dplyr::arrange(
      dplyr::desc(score),
      dplyr::desc(duration),
      dplyr::desc(total_drop)
    )
  
  list(
    run      = scored[1, ],
    slopes   = slopes,
    all_runs = scored
  )
}


find_main_curve_start <- function(time,
                                  o2_fit,
                                  idx_run_start,
                                  pretrans_slope_eps = PRETRANS_SLOPE_EPS,
                                  pretrans_min_run   = PRETRANS_MIN_RUN,
                                  backoff_points     = START_BACKOFF_POINTS,
                                  slope_trend_eps    = SLOPE_TREND_EPS,
                                  smooth_k           = START_SMOOTH_K) {
  slopes <- get_slopes(time, o2_fit)
  slopes_sm <- smooth_slopes(slopes, k = smooth_k)
  
  n <- length(o2_fit)
  
  idx_run_start <- max(3L, min(idx_run_start, n))
  
  idx_seed <- max(2L, min(idx_run_start - backoff_points, idx_run_start))
  left_bound <- max(2L, idx_seed - START_LOOKBACK_POINTS)
  
  search_points <- seq(left_bound, idx_seed - 1L)
  
  best_start <- which.max(o2_fit[left_bound:idx_seed]) + left_bound - 1L
  
  if (length(search_points) == 0 || length(slopes_sm) < pretrans_min_run) {
    return(as.integer(best_start))
  }
  
  n_early <- max(3L, floor(length(slopes_sm) * START_NOISE_FRAC))
  early_idx <- seq_len(n_early)
  
  slope_noise <- stats::mad(slopes_sm[early_idx], na.rm = TRUE)
  
  if (!is.finite(slope_noise) || slope_noise <= 0) {
    slope_noise <- pretrans_slope_eps
  }
  
  real_decline_eps <- max(pretrans_slope_eps, START_NOISE_MULT * slope_noise)
  
  decline_onset <- NA_integer_
  
  for (p in search_points) {
    j1 <- p
    j2 <- min(length(slopes_sm), p + pretrans_min_run - 1L)
    
    if ((j2 - j1 + 1L) < pretrans_min_run) {
      next
    }
    
    seg_s <- slopes_sm[j1:j2]
    
    y1 <- p
    y2 <- min(n, p + pretrans_min_run)
    seg_y <- o2_fit[y1:y2]
    
    if (any(!is.finite(seg_s)) || any(!is.finite(seg_y))) {
      next
    }
    
    cond_neg <- mean(seg_s < -real_decline_eps, na.rm = TRUE) >= START_NEG_FRAC
    
    local_drop <- seg_y[1] - seg_y[length(seg_y)]
    cond_drop <- is.finite(local_drop) && local_drop > (MIN_RANGE_O2 * START_DROP_FRAC)
    
    cond_not_rising <- mean(diff(seg_y) <= START_WIGGLE_TOL, na.rm = TRUE) >=
      START_NOT_RISING_FRAC
    
    if (cond_neg && cond_drop && cond_not_rising) {
      decline_onset <- p
      break
    }
  }
  
  if (is.finite(decline_onset) && !is.na(decline_onset)) {
    peak_left <- max(left_bound, decline_onset - START_LOOKBACK_POINTS)
    peak_right <- min(idx_seed, decline_onset + pretrans_min_run)
    
    peak_window <- peak_left:peak_right
    peak_window <- peak_window[is.finite(o2_fit[peak_window])]
    
    if (length(peak_window) > 0) {
      best_start <- peak_window[which.max(o2_fit[peak_window])]
    } else {
      best_start <- decline_onset
    }
  }
  
  best_start <- max(1L, min(best_start, idx_run_start))
  
  as.integer(best_start)
}


# --------------------------------------------------------------------------
# First-inflection finder.
#
# Returns the index of the FIRST strict local minimum of the smoothed slope
# after idx_peak. The candidate must:
#   - be clearly negative;
#   - occur at least a few points after the start;
#   - occur at least a minimum duration after the start;
#   - be lower than both its left and right local windows.
# --------------------------------------------------------------------------

find_first_inflection <- function(slopes_sm,
                                  time,
                                  idx_peak,
                                  win = END_INFLECTION_WIN,
                                  neg_eps = END_INFLECTION_NEG_EPS,
                                  min_points_after_start = END_MIN_POINTS_AFTER_START,
                                  min_duration_after_start = END_MIN_DURATION_AFTER_START) {
  n <- length(slopes_sm)
  
  if (n < 2 * win + 1) {
    return(NA_integer_)
  }
  
  start_i <- max(
    idx_peak + min_points_after_start,
    idx_peak + win + 1L,
    win + 1L
  )
  
  end_i <- n - win
  
  if (start_i > end_i) {
    return(NA_integer_)
  }
  
  start_time <- time[idx_peak]
  
  for (i in start_i:end_i) {
    s_here <- slopes_sm[i]
    
    if (!is.finite(s_here)) {
      next
    }
    
    if (s_here >= -neg_eps) {
      next
    }
    
    this_time <- time[min(i + 1L, length(time))]
    
    if ((this_time - start_time) < min_duration_after_start) {
      next
    }
    
    left  <- slopes_sm[(i - win):(i - 1L)]
    right <- slopes_sm[(i + 1L):(i + win)]
    
    if (any(!is.finite(c(left, right)))) {
      next
    }
    
    is_local_min <- s_here <= min(left) && s_here <= min(right)
    
    if (is_local_min) {
      return(as.integer(i))
    }
  }
  
  NA_integer_
}


find_end_after_steepest_drop <- function(time,
                                         o2_fit,
                                         idx_peak,
                                         min_gap = END_MIN_GAP,
                                         run_len = END_RECOV_RUN,
                                         smooth_k = END_SMOOTH_K,
                                         recovery_frac = END_RECOV_FRAC,
                                         mono_frac = END_MONO_FRAC,
                                         slope_eps = END_SLOPE_EPS) {
  slopes <- diff(o2_fit) / diff(time)
  slopes[!is.finite(slopes)] <- NA_real_
  
  nsl <- length(slopes)
  
  if (!is.finite(idx_peak) || is.na(idx_peak) || idx_peak >= nsl) {
    return(list(
      idx_end      = NA_integer_,
      idx_steepest = NA_integer_,
      slopes_sm    = slopes,
      end_method   = "invalid_start"
    ))
  }
  
  slopes_sm <- slopes
  
  if (length(slopes) >= smooth_k) {
    sm <- zoo::rollmedian(slopes, k = smooth_k, fill = NA, align = "center")
    sm[is.na(sm)] <- slopes[is.na(sm)]
    slopes_sm <- sm
  }
  
  search_idx <- seq(idx_peak, nsl)
  search_idx <- search_idx[is.finite(slopes_sm[search_idx])]
  
  if (!length(search_idx)) {
    return(list(
      idx_end      = NA_integer_,
      idx_steepest = NA_integer_,
      slopes_sm    = slopes_sm,
      end_method   = "no_valid_slopes"
    ))
  }
  
  idx_steepest <- search_idx[which.min(slopes_sm[search_idx])]
  min_slope <- slopes_sm[idx_steepest]
  
  if (!is.finite(min_slope)) {
    return(list(
      idx_end      = NA_integer_,
      idx_steepest = idx_steepest,
      slopes_sm    = slopes_sm,
      end_method   = "invalid_steepest"
    ))
  }
  
  # Preferred behaviour: end at the first meaningful slope inflection.
  if (isTRUE(END_AT_STEEPEST)) {
    if (isTRUE(END_USE_FIRST_INFLECTION)) {
      idx_first <- find_first_inflection(
        slopes_sm                = slopes_sm,
        time                     = time,
        idx_peak                 = idx_peak,
        win                      = END_INFLECTION_WIN,
        neg_eps                  = END_INFLECTION_NEG_EPS,
        min_points_after_start   = END_MIN_POINTS_AFTER_START,
        min_duration_after_start = END_MIN_DURATION_AFTER_START
      )
      
      if (is.finite(idx_first) && !is.na(idx_first)) {
        return(list(
          idx_end      = as.integer(min(length(o2_fit), idx_first + 1L)),
          idx_steepest = as.integer(idx_first),
          slopes_sm    = slopes_sm,
          end_method   = "first_inflection"
        ))
      }
    }
    
    return(list(
      idx_end      = as.integer(min(length(o2_fit), idx_steepest + 1L)),
      idx_steepest = as.integer(idx_steepest),
      slopes_sm    = slopes_sm,
      end_method   = "global_steepest"
    ))
  }
  
  # Older recovery-after-steepest behaviour.
  target <- min_slope * recovery_frac
  
  start_i <- idx_steepest + min_gap
  end_i   <- nsl - run_len + 1L
  
  if (start_i > end_i) {
    return(list(
      idx_end      = min(length(o2_fit), idx_steepest + 1L),
      idx_steepest = idx_steepest,
      slopes_sm    = slopes_sm,
      end_method   = "fallback_steepest_no_recovery_window"
    ))
  }
  
  for (i in seq(start_i, end_i)) {
    seg <- slopes_sm[i:(i + run_len - 1L)]
    
    if (!all(is.finite(seg))) {
      next
    }
    
    cond_level <- mean(seg >= target, na.rm = TRUE) >= 0.8
    cond_trend <- mean(diff(seg) >= -slope_eps, na.rm = TRUE) >= mono_frac
    
    if (cond_level && cond_trend) {
      return(list(
        idx_end      = as.integer(i + 1L),
        idx_steepest = idx_steepest,
        slopes_sm    = slopes_sm,
        end_method   = "recovery_after_steepest"
      ))
    }
  }
  
  tail_idx <- seq(idx_steepest + 1L, nsl)
  tail_idx <- tail_idx[is.finite(slopes_sm[tail_idx])]
  
  cand <- tail_idx[slopes_sm[tail_idx] >= target]
  
  if (length(cand)) {
    return(list(
      idx_end      = as.integer(cand[1] + 1L),
      idx_steepest = idx_steepest,
      slopes_sm    = slopes_sm,
      end_method   = "first_recovery_candidate"
    ))
  }
  
  list(
    idx_end      = min(length(o2_fit), idx_steepest + 1L),
    idx_steepest = idx_steepest,
    slopes_sm    = slopes_sm,
    end_method   = "fallback_global_steepest"
  )
}


cut_long_plateau_tail <- function(time,
                                  o2_fit,
                                  idx_start,
                                  idx_end,
                                  slope_eps = PLATEAU_SLOPE_EPS,
                                  run_points = PLATEAU_RUN_POINTS,
                                  min_after_start = PLATEAU_MIN_AFTER_START,
                                  keep_points = PLATEAU_KEEP_POINTS) {
  if (isTRUE(END_AT_STEEPEST)) {
    return(as.integer(idx_end))
  }
  
  if (!isTRUE(CUT_LONG_PLATEAU_TAIL)) {
    return(as.integer(idx_end))
  }
  
  n <- length(o2_fit)
  
  idx_start <- max(1L, min(idx_start, n))
  idx_end   <- max(idx_start + 1L, min(idx_end, n))
  
  slopes_abs <- abs(diff(o2_fit) / diff(time))
  slopes_abs[!is.finite(slopes_abs)] <- NA_real_
  
  if (length(slopes_abs) < run_points) {
    return(as.integer(idx_end))
  }
  
  start_i <- max(idx_start + min_after_start, 1L)
  end_i <- min(idx_end - run_points, length(slopes_abs) - run_points + 1L)
  
  if (start_i > end_i) {
    return(as.integer(idx_end))
  }
  
  for (i in seq(start_i, end_i)) {
    seg <- slopes_abs[i:(i + run_points - 1L)]
    
    if (!all(is.finite(seg))) {
      next
    }
    
    cond_plateau <- mean(seg <= slope_eps, na.rm = TRUE) >= 0.80
    
    if (cond_plateau) {
      return(as.integer(min(idx_end, i + keep_points)))
    }
  }
  
  as.integer(idx_end)
}


# --------------------------------------------------------------------------
# Standardised start: argmax of the smoothed spline within
# [START_SEARCH_MIN_TIME, START_SEARCH_MAX_TIME] minutes.
# Falls back to the first finite point if the window is empty.
# --------------------------------------------------------------------------

find_start_simple <- function(time,
                              o2_fit,
                              t_min = START_SEARCH_MIN_TIME,
                              t_max = START_SEARCH_MAX_TIME) {
  n <- length(o2_fit)
  if (n == 0) return(NA_integer_)

  in_win <- which(
    is.finite(o2_fit) &
      is.finite(time) &
      time >= t_min &
      time <= t_max
  )

  if (length(in_win) == 0) {
    first_ok <- which(is.finite(o2_fit) & is.finite(time))
    if (length(first_ok) == 0) return(NA_integer_)
    return(as.integer(first_ok[1]))
  }

  as.integer(in_win[which.max(o2_fit[in_win])])
}


# --------------------------------------------------------------------------
# Standardised end: whichever comes first of
#   (a) idx_start_time + END_MAX_WINDOW_MIN,
#   (b) the slope-recovery point - the first time after the steepest descent
#       that the smoothed slope recovers to recover_frac of its most negative
#       value, i.e. the decline has decelerated and the floor is beginning,
#   (c) the first time the smoothed spline has fallen by END_DROP_FRAC of
#       (O2_fit(start) - min(O2_fit after start))  [safety cap],
#   (d) the end of the recording.
# Returns the index of the chosen end point and a label describing which
# rule fired.
# --------------------------------------------------------------------------

find_end_hybrid <- function(time,
                            o2_fit,
                            idx_start,
                            max_window_min = END_MAX_WINDOW_MIN,
                            drop_frac      = END_DROP_FRAC,
                            recover_frac   = END_SLOPE_RECOVER_FRAC,
                            smooth_k       = END_SMOOTH_K) {
  n <- length(o2_fit)

  if (!is.finite(idx_start) || idx_start < 1L || idx_start >= n) {
    return(list(idx_end = NA_integer_, end_method = "invalid_start"))
  }

  after_idx <- (idx_start + 1L):n
  after_idx <- after_idx[is.finite(o2_fit[after_idx]) & is.finite(time[after_idx])]

  if (length(after_idx) == 0) {
    return(list(idx_end = NA_integer_, end_method = "no_points_after_start"))
  }

  y_start <- o2_fit[idx_start]
  y_min   <- min(o2_fit[after_idx], na.rm = TRUE)
  total_drop <- y_start - y_min

  start_time <- time[idx_start]

  # Rule (a): hard time cap.
  cap_time <- start_time + max_window_min
  idx_cap_candidates <- after_idx[time[after_idx] <= cap_time]
  idx_cap <- if (length(idx_cap_candidates) > 0) {
    max(idx_cap_candidates)
  } else {
    after_idx[1]
  }

  # Rule (b): slope-recovery / flattening onset.
  # Smooth the slope of the spline, find the steepest (most negative) slope
  # after the tip, then take the first point past it where the slope has
  # recovered to recover_frac of that steepest value. That marks where active
  # respiration gives way to the flat (often noisy) floor.
  idx_flat <- NA_integer_
  slopes <- diff(o2_fit) / diff(time)
  slopes[!is.finite(slopes)] <- NA_real_
  if (length(slopes) >= smooth_k) {
    sm <- zoo::rollmedian(slopes, k = smooth_k, fill = NA, align = "center")
    sm[is.na(sm)] <- slopes[is.na(sm)]
    slopes <- sm
  }
  slope_idx <- idx_start:(n - 1L)
  slope_idx <- slope_idx[is.finite(slopes[slope_idx])]
  if (length(slope_idx) > 0) {
    idx_steep <- slope_idx[which.min(slopes[slope_idx])]
    min_slope <- slopes[idx_steep]
    if (is.finite(min_slope) && min_slope < 0) {
      thresh <- recover_frac * min_slope           # less-negative threshold
      after_steep <- (idx_steep + 1L):(n - 1L)
      after_steep <- after_steep[is.finite(slopes[after_steep])]
      rec <- after_steep[slopes[after_steep] >= thresh]
      if (length(rec) > 0) {
        idx_flat <- as.integer(min(rec[1] + 1L, n))
      }
    }
  }

  # Rule (c): drop-fraction threshold (only meaningful if the curve drops).
  idx_drop <- NA_integer_
  if (is.finite(total_drop) && total_drop > 0) {
    target_y <- y_start - drop_frac * total_drop
    hit <- after_idx[o2_fit[after_idx] <= target_y]
    if (length(hit) > 0) {
      idx_drop <- as.integer(hit[1])
    }
  }

  # Rule (d): end of recording.
  idx_eor <- as.integer(n)

  candidates <- c(idx_cap = idx_cap, idx_flat = idx_flat,
                  idx_drop = idx_drop, idx_eor = idx_eor)
  candidates <- candidates[is.finite(candidates)]

  idx_end <- as.integer(min(candidates))

  end_method <- switch(
    names(which.min(candidates))[1],
    idx_cap  = "hybrid_time_cap",
    idx_flat = "hybrid_slope_flatten",
    idx_drop = "hybrid_drop_fraction",
    idx_eor  = "hybrid_end_of_recording",
    "hybrid_unknown"
  )

  list(idx_end = idx_end, end_method = end_method)
}


trim_one_series <- function(df) {
  df <- df %>%
    dplyr::arrange(Time)

  if (nrow(df) < MIN_ROWS_SERIES) {
    return(list(ok = FALSE, reason = "Too few rows", data = NULL, meta = NULL))
  }

  o2_fit <- safe_spline_fit(df$Time, df$Oxygen, spar = TRIM_SPAR)

  if (all(!is.finite(o2_fit))) {
    return(list(ok = FALSE, reason = "Spline failed", data = NULL, meta = NULL))
  }

  if ((max(o2_fit, na.rm = TRUE) - min(o2_fit, na.rm = TRUE)) < MIN_RANGE_O2) {
    return(list(ok = FALSE, reason = "Too flat", data = NULL, meta = NULL))
  }

  df <- df %>%
    dplyr::mutate(O2_fit = o2_fit)

  run_out <- find_main_descending_run(df$Time, df$O2_fit)
  
  main_run <- run_out$run
  slopes   <- run_out$slopes
  all_runs <- run_out$all_runs
  
  if (is.null(main_run) || nrow(main_run) == 0) {
    return(list(ok = FALSE, reason = "No main descending run found", data = NULL, meta = NULL))
  }
  
  idx_run_start <- as.integer(main_run$start[1])
  idx_run_end   <- as.integer(main_run$end[1])

  slopes_sm <- run_out$slopes

  if (isTRUE(END_USE_HYBRID)) {
    # ---- Standardised hybrid path -------------------------------------------
    idx_peak <- find_start_simple(
      time   = df$Time,
      o2_fit = df$O2_fit,
      t_min  = START_SEARCH_MIN_TIME,
      t_max  = START_SEARCH_MAX_TIME
    )

    if (!is.finite(idx_peak) || is.na(idx_peak)) {
      idx_peak <- idx_run_start
    }

    idx_peak_raw <- idx_peak

    end_out <- find_end_hybrid(
      time           = df$Time,
      o2_fit         = df$O2_fit,
      idx_start      = idx_peak,
      max_window_min = END_MAX_WINDOW_MIN,
      drop_frac      = END_DROP_FRAC
    )

    idx_end    <- end_out$idx_end
    end_reason <- end_out$end_method

    if (!is.finite(idx_end) || is.na(idx_end)) {
      idx_end <- min(nrow(df), idx_run_end + 1L)
      end_reason <- "fallback_main_run_end"
    }

    # Approximate "steepest drop" for diagnostic columns only.
    idx_steepest <- {
      slopes_local <- diff(df$O2_fit) / diff(df$Time)
      slopes_local[!is.finite(slopes_local)] <- NA_real_
      after <- seq(idx_peak, length(slopes_local))
      after <- after[is.finite(slopes_local[after])]
      if (length(after)) as.integer(after[which.min(slopes_local[after])]) else NA_integer_
    }
  } else {
    # ---- Legacy inflection / steepest-drop path -----------------------------
    idx_peak <- find_main_curve_start(
      time               = df$Time,
      o2_fit             = df$O2_fit,
      idx_run_start      = idx_run_start,
      pretrans_slope_eps = PRETRANS_SLOPE_EPS,
      pretrans_min_run   = PRETRANS_MIN_RUN,
      backoff_points     = START_BACKOFF_POINTS,
      slope_trend_eps    = SLOPE_TREND_EPS,
      smooth_k           = START_SMOOTH_K
    )

    if (!is.finite(idx_peak) || is.na(idx_peak)) {
      idx_peak <- idx_run_start
    }

    idx_peak_raw <- idx_peak

    end_out <- find_end_after_steepest_drop(
      time          = df$Time,
      o2_fit        = df$O2_fit,
      idx_peak      = idx_peak,
      min_gap       = END_MIN_GAP,
      run_len       = END_RECOV_RUN,
      smooth_k      = END_SMOOTH_K,
      recovery_frac = END_RECOV_FRAC,
      mono_frac     = END_MONO_FRAC,
      slope_eps     = END_SLOPE_EPS
    )

    idx_end      <- end_out$idx_end
    idx_steepest <- end_out$idx_steepest
    slopes_sm    <- end_out$slopes_sm
    end_method   <- end_out$end_method

    if (!is.finite(idx_end) || is.na(idx_end)) {
      idx_end <- min(nrow(df), idx_run_end + 1L)
      end_reason <- "fallback_main_run_end"
    } else {
      end_reason <- end_method
    }

    idx_end_before_plateau_cut <- idx_end

    idx_end <- cut_long_plateau_tail(
      time            = df$Time,
      o2_fit          = df$O2_fit,
      idx_start       = idx_peak,
      idx_end         = idx_end,
      slope_eps       = PLATEAU_SLOPE_EPS,
      run_points      = PLATEAU_RUN_POINTS,
      min_after_start = PLATEAU_MIN_AFTER_START,
      keep_points     = PLATEAU_KEEP_POINTS
    )

    if (idx_end < idx_end_before_plateau_cut) {
      end_reason <- paste0(end_reason, "+plateau_tail_cut")
    }
  }

  idx_peak <- max(1L, min(idx_peak, nrow(df)))
  idx_end  <- max(idx_peak + 1L, min(idx_end, nrow(df)))

  if (!is.finite(idx_steepest) || is.na(idx_steepest)) {
    idx_steepest <- max(idx_peak, min(nrow(df) - 1L, idx_run_start))
  }

  if (idx_end <= idx_peak) {
    return(list(ok = FALSE, reason = "Invalid peak/end order", data = NULL, meta = NULL))
  }

  # ---- Hybrid guards: minimum kept duration and point count -----------------
  if (isTRUE(END_USE_HYBRID)) {
    n_kept    <- idx_end - idx_peak + 1L
    dur_kept  <- df$Time[idx_end] - df$Time[idx_peak]

    if (n_kept < MIN_FIT_POINTS) {
      return(list(
        ok = FALSE,
        reason = sprintf("Kept window has %d points (< MIN_FIT_POINTS = %d)",
                         n_kept, MIN_FIT_POINTS),
        data = NULL, meta = NULL
      ))
    }

    if (!is.finite(dur_kept) || dur_kept < MIN_FIT_DURATION_MIN) {
      return(list(
        ok = FALSE,
        reason = sprintf("Kept window is %.2f min (< MIN_FIT_DURATION_MIN = %.2f)",
                         dur_kept, MIN_FIT_DURATION_MIN),
        data = NULL, meta = NULL
      ))
    }
  }
  
  trimmed <- df[idx_peak:idx_end, ] %>%
    dplyr::mutate(
      curve_code           = df$curve_code[1],
      series_id            = df$series_id[1],
      peak_time_raw        = df$Time[idx_peak_raw],
      peak_time            = df$Time[idx_peak],
      main_run_start_time  = df$Time[idx_run_start],
      main_run_end_time    = df$Time[min(idx_run_end + 1L, nrow(df))],
      steepest_drop_time   = df$Time[min(idx_steepest + 1L, nrow(df))],
      end_time_chosen      = df$Time[idx_end],
      peak_idx_raw         = idx_peak_raw,
      peak_idx             = idx_peak,
      main_run_start_idx   = idx_run_start,
      main_run_end_idx     = idx_run_end,
      steepest_drop_idx    = idx_steepest,
      end_idx              = idx_end,
      start_shift_points   = idx_run_start - idx_peak,
      start_shift_time     = df$Time[idx_run_start] - df$Time[idx_peak],
      end_reason           = end_reason
    )
  
  meta <- tibble::tibble(
    curve_code = df$curve_code[1],
    series_id = df$series_id[1],
    T = df$T[1],
    OTU = df$OTU[1],
    Replicate = df$Replicate[1],
    n_points_total = nrow(df),
    n_points_trimmed = nrow(trimmed),
    peak_time_raw = df$Time[idx_peak_raw],
    peak_time = df$Time[idx_peak],
    main_run_start_time = df$Time[idx_run_start],
    main_run_end_time = df$Time[min(idx_run_end + 1L, nrow(df))],
    steepest_drop_time = df$Time[min(idx_steepest + 1L, nrow(df))],
    chosen_end_time = df$Time[idx_end],
    fit_duration_min = df$Time[idx_end] - df$Time[idx_peak],
    delta_Ninoc_to_N0_min = df$Time[idx_peak],
    start_shift_points = idx_run_start - idx_peak,
    start_shift_time = df$Time[idx_run_start] - df$Time[idx_peak],
    end_reason = end_reason,
    main_run_duration = main_run$duration[1],
    main_run_drop = main_run$total_drop[1],
    main_run_score = main_run$score[1]
  )
  
  list(
    ok            = TRUE,
    reason        = NA_character_,
    data          = trimmed,
    meta          = meta,
    idx_peak_raw  = idx_peak_raw,
    idx_peak      = idx_peak,
    idx_run_start = idx_run_start,
    idx_run_end   = idx_run_end,
    idx_steepest  = idx_steepest,
    idx_end       = idx_end,
    slopes        = slopes,
    slopes_sm     = slopes_sm,
    all_runs      = all_runs,
    full_df       = df,
    manual_override = FALSE,
    manual_keep_start = NA_real_,
    manual_keep_end   = NA_real_
  )
}


# =============================================================================
# 5) Load data
# =============================================================================

if (!file.exists(LONG_CSV)) {
  stop("Input file not found: ", LONG_CSV)
}

raw_before_cut <- readr::read_csv(LONG_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    T         = as.numeric(T),
    OTU     = as.integer(OTU),
    Replicate = as.character(Replicate),
    Time      = as.numeric(Time),
    Oxygen    = as.numeric(Oxygen),
    series_id = paste0("T=", T, " | OTU=", OTU, " | Rep=", Replicate)
  )

required_cols <- c("Time", "T", "OTU", "Replicate", "Oxygen")
missing_cols <- setdiff(required_cols, names(raw_before_cut))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in Oxygen_All_Long.csv: ",
    paste(missing_cols, collapse = ", ")
  )
}

raw <- raw_before_cut %>%
  dplyr::filter(Time <= MAX_TIME_MIN)

if (nrow(raw) == 0) {
  stop("No rows remain after applying MAX_TIME_MIN = ", MAX_TIME_MIN, ".")
}

message(
  "Time cutoff applied: kept ", nrow(raw), " of ", nrow(raw_before_cut),
  " rows with Time <= ", MAX_TIME_MIN, " min. Max retained Time = ",
  max(raw$Time, na.rm = TRUE)
)


# =============================================================================
# 6) Create short code for each curve
# =============================================================================

curve_key <- raw %>%
  dplyr::distinct(series_id, T, OTU, Replicate) %>%
  dplyr::mutate(
    rep_num = suppressWarnings(as.integer(stringr::str_remove(Replicate, "^R")))
  ) %>%
  dplyr::arrange(T, OTU, rep_num, Replicate) %>%
  dplyr::mutate(
    curve_code = paste0(
      "C",
      stringr::str_pad(dplyr::row_number(), width = 3, pad = "0")
    )
  ) %>%
  dplyr::select(curve_code, series_id, T, OTU, Replicate)

raw <- raw %>%
  dplyr::left_join(
    curve_key,
    by = c("series_id", "T", "OTU", "Replicate")
  )

readr::write_csv(
  curve_key,
  file.path(tables_dir, "Oxygen_Curve_Code_Key.csv")
)


# =============================================================================
# 7) Validate manual overrides
# =============================================================================

if (nrow(manual_trim_overrides) > 0) {
  needed_override_cols <- c("curve_code", "keep_start", "keep_end")
  missing_override_cols <- setdiff(needed_override_cols, names(manual_trim_overrides))
  
  if (length(missing_override_cols) > 0) {
    stop(
      "manual_trim_overrides is missing columns: ",
      paste(missing_override_cols, collapse = ", ")
    )
  }
  
  bad_codes <- setdiff(manual_trim_overrides$curve_code, curve_key$curve_code)
  
  if (length(bad_codes) > 0) {
    stop(
      "These curve_code values were not found: ",
      paste(bad_codes, collapse = ", ")
    )
  }
  
  bad_ranges <- manual_trim_overrides %>%
    dplyr::filter(
      (!is.na(keep_start) & !is.finite(keep_start)) |
        (!is.na(keep_end) & !is.finite(keep_end))
    )
  
  if (nrow(bad_ranges) > 0) {
    stop("Some manual overrides have invalid keep_start or keep_end values.")
  }
}


# =============================================================================
# 8) Trim each series
# =============================================================================

series_ids <- unique(raw$series_id)

trimmed_lst <- vector("list", length(series_ids))
names(trimmed_lst) <- series_ids

meta_lst <- vector("list", length(series_ids))
names(meta_lst) <- series_ids

diag_lst <- vector("list", length(series_ids))
names(diag_lst) <- series_ids

skipped_log <- tibble::tibble(
  curve_code = character(),
  series_id  = character(),
  reason     = character()
)

for (sid in series_ids) {
  df <- raw %>%
    dplyr::filter(series_id == sid) %>%
    dplyr::arrange(Time)
  
  this_code <- df$curve_code[1]
  
  out <- trim_one_series(df)
  
  if (!isTRUE(out$ok)) {
    skipped_log <- skipped_log %>%
      tibble::add_row(
        curve_code = this_code,
        series_id = sid,
        reason = out$reason
      )
    
    next
  }
  
  ov <- manual_trim_overrides %>%
    dplyr::filter(curve_code == this_code)
  
  if (nrow(ov) > 0) {
    auto_start <- df$Time[out$idx_peak]
    auto_end   <- df$Time[out$idx_end]
    
    final_start <- if (!is.na(ov$keep_start[1])) ov$keep_start[1] else auto_start
    final_end   <- if (!is.na(ov$keep_end[1]))   ov$keep_end[1]   else auto_end
    
    final_end <- min(final_end, max(df$Time, na.rm = TRUE), MAX_TIME_MIN, na.rm = TRUE)
    
    if (!is.finite(final_start) || !is.finite(final_end) || final_end <= final_start) {
      skipped_log <- skipped_log %>%
        tibble::add_row(
          curve_code = this_code,
          series_id = sid,
          reason = "Manual override gave invalid final start/end"
        )
      
      next
    }
    
    manual_trimmed <- df %>%
      dplyr::filter(Time >= final_start, Time <= final_end)
    
    if (nrow(manual_trimmed) < 2) {
      skipped_log <- skipped_log %>%
        tibble::add_row(
          curve_code = this_code,
          series_id = sid,
          reason = "Manual override produced fewer than 2 rows"
        )
      
      next
    }
    
    manual_trimmed <- manual_trimmed %>%
      dplyr::mutate(
        O2_fit = safe_spline_fit(Time, Oxygen, spar = TRIM_SPAR),
        curve_code = this_code,
        series_id = df$series_id[1],
        peak_time_raw = auto_start,
        peak_time = final_start,
        main_run_start_time = df$Time[out$idx_run_start],
        main_run_end_time = df$Time[min(out$idx_run_end + 1L, nrow(df))],
        steepest_drop_time = if (is.finite(out$idx_steepest)) {
          df$Time[min(out$idx_steepest + 1L, nrow(df))]
        } else {
          NA_real_
        },
        end_time_chosen = final_end,
        peak_idx_raw = out$idx_peak_raw,
        peak_idx = NA_integer_,
        main_run_start_idx = out$idx_run_start,
        main_run_end_idx = out$idx_run_end,
        steepest_drop_idx = out$idx_steepest,
        end_idx = NA_integer_,
        start_shift_points = NA_integer_,
        start_shift_time = final_start - auto_start,
        end_reason = "manual_override_partial_or_full"
      )
    
    manual_meta <- tibble::tibble(
      curve_code = this_code,
      series_id = df$series_id[1],
      T = df$T[1],
      OTU = df$OTU[1],
      Replicate = df$Replicate[1],
      n_points_total = nrow(df),
      n_points_trimmed = nrow(manual_trimmed),
      peak_time_raw = auto_start,
      peak_time = final_start,
      main_run_start_time = df$Time[out$idx_run_start],
      main_run_end_time = df$Time[min(out$idx_run_end + 1L, nrow(df))],
      steepest_drop_time = if (is.finite(out$idx_steepest)) {
        df$Time[min(out$idx_steepest + 1L, nrow(df))]
      } else {
        NA_real_
      },
      chosen_end_time = final_end,
      fit_duration_min = final_end - final_start,
      delta_Ninoc_to_N0_min = final_start,
      start_shift_points = NA_integer_,
      start_shift_time = final_start - auto_start,
      end_reason = "manual_override_partial_or_full",
      main_run_duration = if (!is.null(out$all_runs) && nrow(out$all_runs) > 0) {
        out$all_runs$duration[1]
      } else {
        NA_real_
      },
      main_run_drop = if (!is.null(out$all_runs) && nrow(out$all_runs) > 0) {
        out$all_runs$total_drop[1]
      } else {
        NA_real_
      },
      main_run_score = if (!is.null(out$all_runs) && nrow(out$all_runs) > 0) {
        out$all_runs$score[1]
      } else {
        NA_real_
      }
    )
    
    out$data <- manual_trimmed
    out$meta <- manual_meta
    out$manual_override <- TRUE
    out$manual_keep_start <- if (!is.na(ov$keep_start[1])) ov$keep_start[1] else NA_real_
    out$manual_keep_end   <- if (!is.na(ov$keep_end[1])) ov$keep_end[1] else NA_real_
    out$final_start <- final_start
    out$final_end   <- final_end
  } else {
    out$final_start <- df$Time[out$idx_peak]
    out$final_end   <- min(df$Time[out$idx_end], MAX_TIME_MIN)
  }
  
  trimmed_lst[[sid]] <- out$data
  meta_lst[[sid]]    <- out$meta
  diag_lst[[sid]]    <- out
}


# =============================================================================
# 9) Save trimmed data and metadata
# =============================================================================

trimmed <- dplyr::bind_rows(trimmed_lst)
trim_meta <- dplyr::bind_rows(meta_lst)

if (nrow(trimmed) == 0) {
  stop(
    "No curves were successfully trimmed. Check Skipped_Series_Log.csv and trimming parameters."
  )
}

trimmed <- trimmed %>%
  dplyr::filter(Time <= MAX_TIME_MIN)

if (max(trimmed$Time, na.rm = TRUE) > MAX_TIME_MIN) {
  stop("Internal error: trimmed output still contains Time > MAX_TIME_MIN.")
}

readr::write_csv(
  trimmed,
  file.path(tables_dir, "Oxygen_Data_Smoothed_Trimmed.csv")
)

filtered <- trimmed %>%
  dplyr::select(curve_code, T, OTU, Replicate, Time, Oxygen, O2_fit)

readr::write_csv(filtered, IN_CSV)

readr::write_csv(
  skipped_log,
  file.path(tables_dir, "Skipped_Series_Log.csv")
)

trim_meta <- trim_meta %>%
  dplyr::mutate(
    rep_num = suppressWarnings(as.integer(stringr::str_remove(Replicate, "^R")))
  ) %>%
  dplyr::arrange(T, OTU, rep_num, Replicate) %>%
  dplyr::select(-rep_num)

readr::write_csv(trim_meta, TRIM_META_CSV)


# =============================================================================
# 10) Diagnostics PDF
# =============================================================================

pdf_file <- file.path(figures_dir, "oxygen_trimming_diagnostics.pdf")

if (file.exists(pdf_file)) {
  file.remove(pdf_file)
}

pdf(
  pdf_file,
  width = 8.2,
  height = 6.2
)

for (sid in names(diag_lst)) {
  out <- diag_lst[[sid]]
  
  if (is.null(out)) {
    next
  }
  
  df <- out$full_df %>%
    dplyr::filter(Time <= MAX_TIME_MIN)
  
  this_code <- df$curve_code[1]
  
  idx_peak_raw  <- out$idx_peak_raw
  idx_peak      <- out$idx_peak
  idx_run_start <- out$idx_run_start
  idx_steepest  <- out$idx_steepest
  idx_end       <- out$idx_end
  
  xmin_show <- max(min(out$final_start, MAX_TIME_MIN), min(df$Time, na.rm = TRUE))
  xmax_show <- min(out$final_end, MAX_TIME_MIN, max(df$Time, na.rm = TRUE))

  # Zoom the panel to the kept window plus a margin, so the long flat tail past
  # the chosen end is cropped out of the display (data are unchanged).
  kept_span  <- max(xmax_show - xmin_show, 0)
  plot_pad   <- max(PLOT_PAD_MIN_MIN, PLOT_PAD_FRAC * kept_span)
  x_plot_lo  <- max(min(df$Time, na.rm = TRUE), xmin_show - plot_pad)
  x_plot_hi  <- min(MAX_TIME_MIN, max(df$Time, na.rm = TRUE), xmax_show + plot_pad)

  rect_df <- tibble::tibble(
    xmin = xmin_show,
    xmax = xmax_show,
    ymin = -Inf,
    ymax = Inf
  )
  
  subtitle_text <- if (isTRUE(out$manual_override)) {
    paste0(
      "MANUAL OVERRIDE | final kept region = ",
      xmin_show, " to ", xmax_show, " min | ",
      "entered start = ",
      ifelse(is.na(out$manual_keep_start), "auto", as.character(out$manual_keep_start)),
      " | entered end = ",
      ifelse(is.na(out$manual_keep_end), "auto", as.character(out$manual_keep_end)),
      " | max time = ", MAX_TIME_MIN,
      " | spar = ", TRIM_SPAR
    )
  } else {
    paste0(
      "Grey dotted = detected oxygen tip | ",
      "Black = chosen start | ",
      "Green dashed = main descending run start | ",
      "Purple dashed = selected inflection/steepest drop | ",
      "Magenta = chosen end | ",
      "End reason = ", out$meta$end_reason[1], " | ",
      "max time = ", MAX_TIME_MIN,
      " | spar = ", TRIM_SPAR
    )
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(Time, Oxygen)) +
    ggplot2::geom_rect(
      data = rect_df,
      inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "orange",
      alpha = 0.12
    ) +
    ggplot2::geom_line(colour = "grey60", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = O2_fit), colour = "blue", linewidth = 1.1) +
    ggplot2::coord_cartesian(xlim = c(x_plot_lo, x_plot_hi))
  
  if (!isTRUE(out$manual_override)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = df$Time[pmin(idx_peak_raw, nrow(df))],
        colour = "grey20",
        linetype = "dotted",
        linewidth = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = df$Time[pmin(idx_peak, nrow(df))],
        colour = "black",
        linetype = "solid",
        linewidth = 0.9
      ) +
      ggplot2::geom_vline(
        xintercept = df$Time[pmin(idx_run_start, nrow(df))],
        colour = "darkgreen",
        linetype = "dashed",
        linewidth = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = df$Time[pmin(idx_steepest + 1L, nrow(df))],
        colour = "purple",
        linetype = "dashed",
        linewidth = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = df$Time[pmin(idx_end, nrow(df))],
        colour = "magenta",
        linetype = "solid",
        linewidth = 0.9
      )
  } else {
    p <- p +
      ggplot2::geom_vline(
        xintercept = xmin_show,
        colour = "black",
        linetype = "solid",
        linewidth = 0.9
      ) +
      ggplot2::geom_vline(
        xintercept = xmax_show,
        colour = "magenta",
        linetype = "solid",
        linewidth = 0.9
      )
  }
  
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = X_MAJOR_BREAKS_N),
      minor_breaks = scales::pretty_breaks(n = X_MAJOR_BREAKS_N * 2)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = Y_MAJOR_BREAKS_N),
      minor_breaks = scales::pretty_breaks(n = Y_MAJOR_BREAKS_N * 2)
    ) +
    ggplot2::labs(
      title = paste0(this_code, " | ", sid),
      subtitle = subtitle_text,
      x = "Time (min)",
      y = "O2 (mg L^-1)"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.ticks.length = grid::unit(0.18, "cm"),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_line(colour = "grey94", linewidth = 0.2)
    )
  
  print(p)
}

dev.off()


# =============================================================================
# 11) Completion messages
# =============================================================================

message("Saved: ", file.path(tables_dir, "Oxygen_Curve_Code_Key.csv"))
message("Saved: ", file.path(tables_dir, "Oxygen_Data_Smoothed_Trimmed.csv"))
message("Saved: ", IN_CSV)
message("Saved: ", TRIM_META_CSV)
message("Saved: ", file.path(tables_dir, "Skipped_Series_Log.csv"))
message("Saved: ", pdf_file)
message(
  "Done: oxygen trimming completed successfully. Max retained time = ",
  max(trimmed$Time, na.rm = TRUE), " min."
)