# ============================================================
# R/revision_2026/empirical/01_data.R
#
# Data access for the empirical application (scripts 02-04).
#
# The HELIUS preprocessing and the perspective paper's estimates live in the
# perspective repository (slow-fast-perspective). Point to it with the
# environment variable SLOWFAST_PERSPECTIVE_DIR; the default assumes the two
# repositories sit side by side:
#
#   parent/
#     slow-fast-computational-model/   <- run scripts from here
#     slow-fast-perspective/           <- R/00-03, data/helius.sav, results/
#
# Two modes, set with the environment variable SLOWFAST_DATA_MODE:
#
#   "helius"     (default) runs the perspective repo's preprocessing scripts
#                00-03 and uses the real HELIUS baseline data.
#
#   "synthetic"  builds a HELIUS-like data set from the perspective paper's
#                published group-level estimates, so that the full pipeline
#                can be run and checked without access to HELIUS. Results in
#                this mode are only for checking the pipeline -- not real results.
#
# Both modes return `sf_data`, a data frame with:
#   symptoms (0/1)  anh dep slp ene app glt con mot sui
#   srl_z           standardized Slow Risk Load
#   high_srl        1 = High Slow Risk (median split), 0 = Low
#   age_z, female, non_dutch   covariates for adjusted models
#   older           1 = older half (median split of age)
#
# Outputs go to res/revision_2026/empirical/<mode>/ in this repository and
# contain only aggregate estimates. The cached synthetic data set is written
# to res/revision_2026/empirical/private/, which is git-ignored.
# ============================================================

sf_symptoms <- c("anh", "dep", "slp", "ene", "app", "glt", "con", "mot", "sui")
sf_find_perspective <- function() {
  env <- Sys.getenv("SLOWFAST_PERSPECTIVE_DIR", "")
  candidates <- if (nzchar(env)) env else
    c("inputs/perspective", "../slow-fast-perspective", "../slow-fast-perspective/slow-fast-perspective")
  for (d in candidates) {
    if (dir.exists(file.path(d, "results", "network_invariance_exact"))) {
      return(normalizePath(d))
    }
  }
  normalizePath(candidates[1], mustWork = FALSE)
}
SF_PERSPECTIVE <- sf_find_perspective()
SF_OUTPUT <- file.path("res", "revision_2026", "empirical")
SF_PRIVATE <- file.path(SF_OUTPUT, "private")
SF_PERSP <- file.path(SF_PERSPECTIVE, "results", "network_invariance_exact", "slow_risk")
if (!dir.exists(SF_PERSP)) {
  stop("Perspective-paper results not found at ", SF_PERSP, ".\n",
       "Set SLOWFAST_PERSPECTIVE_DIR to the slow-fast-perspective repository.")
}
dir.create(SF_PRIVATE, recursive = TRUE, showWarnings = FALSE)

sf_mode <- function() {
  m <- Sys.getenv("SLOWFAST_DATA_MODE", "helius")
  match.arg(m, c("helius", "synthetic"))
}


# ------------------------------------------------------------------------
# Perspective-paper estimates (already in this repository)
# ------------------------------------------------------------------------

sf_read_perspective <- function(model = "equal_J_free_h") {
  d <- file.path(SF_PERSP, "exact_multigroup_ising")
  a <- utils::read.csv(file.path(d, paste0("activation_parameters_", model, ".csv")))
  J1 <- as.matrix(utils::read.csv(file.path(d, paste0("J_", model, "_group_1.csv")), row.names = 1))
  J2 <- as.matrix(utils::read.csv(file.path(d, paste0("J_", model, "_group_2.csv")), row.names = 1))
  list(
    h_low = stats::setNames(a$h[a$group == "Low Slow Risk"], a$symptom[a$group == "Low Slow Risk"]),
    h_high = stats::setNames(a$h[a$group == "High Slow Risk"], a$symptom[a$group == "High Slow Risk"]),
    J_low = J1[sf_symptoms, sf_symptoms],
    J_high = J2[sf_symptoms, sf_symptoms]
  )
}

sf_perspective_observed <- function() {
  d <- file.path(SF_PERSP, "exact_multigroup_ising")
  mc <- utils::read.csv(file.path(d, "model_comparison_exact.csv"))
  ll <- stats::setNames(mc$logLik, mc$model)
  nct <- utils::read.csv(file.path(SF_PERSP, "nct_10000", "nct_summary_10000.csv"))
  dfree <- utils::read.csv(file.path(d, "activation_differences_free_J_free_h.csv"))
  deq <- utils::read.csv(file.path(d, "activation_differences_equal_J_free_h.csv"))
  free <- sf_read_perspective("free_J_free_h")
  ut <- upper.tri(free$J_low)
  data.frame(
    LR_J = 2 * (ll[["free_J_free_h"]] - ll[["equal_J_free_h"]]),
    LR_h = 2 * (ll[["free_J_free_h"]] - ll[["free_J_equal_h"]]),
    bic_best = mc$model[which.min(mc$BIC)],
    aic_best = mc$model[which.min(mc$AIC)],
    dBIC_commonJ_vs_free = mc$BIC[mc$model == "equal_J_free_h"] - mc$BIC[mc$model == "free_J_free_h"],
    dBIC_commonJ_vs_commonh = mc$BIC[mc$model == "equal_J_free_h"] - mc$BIC[mc$model == "free_J_equal_h"],
    sum_dJ_free = sum((free$J_high - free$J_low)[ut]),
    n_dh_free_positive = sum(dfree$h_difference_group2_minus_group1 > 0),
    n_dh_commonJ_positive = sum(deq$h_difference_group2_minus_group1 > 0),
    dh_slp_free = dfree$h_difference_group2_minus_group1[dfree$symptom == "slp"],
    dh_slp_commonJ = deq$h_difference_group2_minus_group1[deq$symptom == "slp"],
    nct_M = nct$network_structure_statistic_M,
    nct_M_p = nct$network_structure_p,
    nct_S = nct$global_strength_difference_S,
    nct_S_p = nct$global_strength_p,
    nct_strength_low = nct$global_strength_group_1,
    nct_strength_high = nct$global_strength_group_2,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------
# Real HELIUS data
# ------------------------------------------------------------------------

sf_load_helius <- function() {
  need <- c("R/00_setup.R", "R/01_prepare_analysis_data.R",
            "R/02_model_helpers.R", "R/03_prepare_binary_data.R")
  if (!all(file.exists(file.path(SF_PERSPECTIVE, need))) ||
      !file.exists(file.path(SF_PERSPECTIVE, "data", "helius.sav"))) {
    stop("HELIUS mode needs ", SF_PERSPECTIVE, "/R/00-03 and data/helius.sav. ",
         "Set SLOWFAST_PERSPECTIVE_DIR, or SLOWFAST_DATA_MODE=synthetic.")
  }
  # The perspective scripts use paths relative to their own repository root.
  old_wd <- setwd(SF_PERSPECTIVE)
  on.exit(setwd(old_wd), add = TRUE)
  env <- new.env()
  for (f in need) sys.source(f, envir = env)
  d <- env$analysis_df_bin
  out <- data.frame(
    d[, sf_symptoms],
    srl_z = d$slow_risk_z,
    high_srl = as.integer(d$slow_group == "High Slow Risk"),
    age = d$age,
    female = as.integer(d$gender_grp == "Female"),
    non_dutch = as.integer(d$dutch_grp == "Non-Dutch"),
    older = as.integer(d$age_grp == "Older")
  )
  # Same analytic sample as the perspective paper: complete symptoms and SRL.
  # Ethnicity can be missing; those rows are dropped only from the
  # ethnicity-adjusted models and the Dutch vs Non-Dutch comparison.
  core <- c(sf_symptoms, "srl_z", "high_srl", "age", "female", "older")
  out <- out[stats::complete.cases(out[, core]), , drop = FALSE]
  out$age_z <- as.numeric(scale(out$age))
  for (v in sf_symptoms) out[[v]] <- as.integer(out[[v]])
  if (any(is.na(out$non_dutch))) {
    message(sum(is.na(out$non_dutch)), " participant(s) with missing ethnicity: kept in ",
            "unadjusted analyses, excluded from adjusted models and the ethnicity comparison.")
  }
  out
}


# ------------------------------------------------------------------------
# Synthetic HELIUS-like data (for checking the pipeline only)
# ------------------------------------------------------------------------
#
# Generating model = the slow-fast cross-sectional model with shared J:
#   h_n = h0 + gamma * srl_z_n
# SRL is right-skewed (longer upper tail, as in the perspective paper's Fig. 3).
# h0, gamma and J are calibrated so that the perspective paper's
# equal_J_free_h analysis of the synthetic data reproduces its published
# Low/High activations and interactions (see sf_calibrate_context()).

sf_skewed_srl <- function(n, shape = 6) {
  g <- stats::rgamma(n, shape = shape)
  as.numeric(scale(g))
}

sf_make_synthetic <- function(n = 23283L, seed = 20261019L) {
  source("R/revision_2026/empirical/00_context_ising_core.R", local = TRUE)
  set.seed(seed)
  st <- ci_setup(sf_symptoms)
  srl_z <- sf_skewed_srl(n)
  cal <- sf_calibrate_context(srl_z, st)
  Z <- cbind(1, srl_z)
  X <- ci_sample(cal$B, Z, st)
  colnames(X) <- sf_symptoms
  age <- pmin(70, pmax(18, stats::rnorm(n, 43.8 + 2 * srl_z, 13.4)))
  non_dutch <- stats::rbinom(n, 1, stats::plogis(1.4 + 1.2 * srl_z))
  out <- data.frame(
    X,
    srl_z = srl_z,
    high_srl = as.integer(srl_z > stats::median(srl_z)),
    age = age,
    female = stats::rbinom(n, 1, .574),
    non_dutch = non_dutch,
    older = as.integer(age >= stats::median(age))
  )
  out$age_z <- as.numeric(scale(out$age))
  attr(out, "generating_B") <- cal$B
  out
}


# ------------------------------------------------------------------------
# Calibration of the shared-J context model to group-level targets
# ------------------------------------------------------------------------
#
# Finds h0, gamma (per symptom) and a shared J such that, given the SRL
# values srl_z, fitting the perspective paper's common-interaction model
# (equal_J_free_h) to the model's *expected* sufficient statistics returns the
# target h_low, h_high and J. Simple fixed-point iteration: the mixture over
# SRL within each group only slightly distorts the fitted parameters.

sf_calibrate_context <- function(srl_z, setup, targets = sf_read_perspective("equal_J_free_h"),
                                 iters = 25L, tol = 1e-4, verbose = FALSE) {
  st <- setup
  hi <- srl_z > stats::median(srl_z)
  mz <- c(mean(srl_z[!hi]), mean(srl_z[hi]))
  tJ <- targets$J_low[st$pair_index]
  # first-order start
  gamma <- (targets$h_high - targets$h_low) / diff(mz)
  h0 <- targets$h_low - gamma * mz[1]
  J <- tJ
  for (it in seq_len(iters)) {
    B <- cbind(c(h0, J), c(gamma, rep(0, st$q)))
    s_lo <- ci_expected_stats(B, cbind(1, srl_z[!hi]), st)
    s_hi <- ci_expected_stats(B, cbind(1, srl_z[hi]), st)
    f <- ci_two_group_suffstat_fit(list(s_lo, s_hi), list(sum(!hi), sum(hi)), st)$equal_J_free_h
    e_lo <- targets$h_low - f$h1; e_hi <- targets$h_high - f$h2; e_J <- tJ - f$J1
    err <- max(abs(c(e_lo, e_hi, e_J)))
    if (verbose) message(sprintf("calibration iter %d  max error %.2e", it, err))
    if (err < tol) break
    # h_g ~ h0 + gamma * mz[g]: solve the 2x2 system for the corrections
    dg <- (e_hi - e_lo) / diff(mz)
    gamma <- gamma + dg
    h0 <- h0 + e_lo - dg * mz[1]
    J <- J + e_J
  }
  list(B = cbind(c(h0, J), c(gamma, rep(0, st$q))),
       h0 = stats::setNames(h0, st$symptoms),
       gamma = stats::setNames(gamma, st$symptoms),
       J = J, max_error = err, iterations = it)
}


# ------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------

sf_load_data <- function() {
  if (sf_mode() == "helius") {
    d <- sf_load_helius()
  } else {
    cache <- file.path(SF_PRIVATE, "synthetic_data.rds")
    if (file.exists(cache)) {
      d <- readRDS(cache)
    } else {
      d <- sf_make_synthetic()
      saveRDS(d, cache)
    }
    message("SLOWFAST_DATA_MODE=synthetic -- fake data, only for checking the pipeline.")
  }
  d
}

sf_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE)
  invisible(file)
}

sf_out <- function(...) {
  # mode-specific output folder, so synthetic tests never overwrite real results
  file.path(SF_OUTPUT, sf_mode(), ...)
}
