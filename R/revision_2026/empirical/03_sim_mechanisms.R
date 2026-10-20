# ============================================================
# R/revision_2026/empirical/03_sim_mechanisms.R
#
# Step 2: Which data-generating mechanism reproduces the perspective paper's
# empirical pattern?
#
# We generate HELIUS-sized cross-sectional samples (one independent draw per
# person, each from that person's exact 512-state distribution) under several
# generating conditions, then run the perspective paper's analysis on each:
# the four exact multigroup Ising models on the Low/High median split, plus
# the network-comparison statistics.
#
# Generating conditions:
#
#   null             one set of h and J for everyone
#   fixed_field      groups differ only in a constant activation shift, shared J
#   context          shared J; each person's activations shift with their own
#                    SRL (slow-fast model). The median split therefore leaves
#                    unmodelled within-group context variation.
#   true_coupling    groups have genuinely different J (and h), no within-group
#                    variation
#   context_coupling SRL shifts activations AND moderates J (step 1's
#                    ctx_coupling fit; only when step 1 has been run)
#
# Parameters come from step 1 when its results exist for the current data
# mode; otherwise the context condition is calibrated to the perspective-paper
# paper's published group-level estimates (sf_calibrate_context). The
# null / fixed_field / true_coupling conditions use the perspective paper's
# fitted equal_J_equal_h / equal_J_free_h / free_J_free_h parameters.
#
# For each condition we report:
#   (a) the population-level result: the four models fitted to the condition's
#       expected sufficient statistics, giving the expected likelihood-ratio
#       non-centrality (no simulation noise);
#   (b) the distribution over replications of the same summaries, and where
#       the observed data fall in that distribution;
#   (c) IsingFit-based NCT statistics (M = max edge difference, S = global
#       strength difference) over a subset of replications, benchmarked
#       against the null condition. Optionally, the full permutation NCT on a
#       few replications.
#
# Run from this repository's root (slow-fast-computational-model):
#   Rscript R/revision_2026/empirical/03_sim_mechanisms.R
#   SLOWFAST_DATA_MODE=synthetic Rscript R/revision_2026/empirical/03_sim_mechanisms.R
# ============================================================

source("R/revision_2026/empirical/00_context_ising_core.R")
source("R/revision_2026/empirical/01_data.R")

N_REPS <- as.integer(Sys.getenv("SLOWFAST_N_REPS", "200"))
N_NETWORK_REPS <- as.integer(Sys.getenv("SLOWFAST_N_NETWORK_REPS", "50"))
N_FULL_NCT_REPS <- as.integer(Sys.getenv("SLOWFAST_N_FULL_NCT_REPS", "0"))
FULL_NCT_IT <- as.integer(Sys.getenv("SLOWFAST_FULL_NCT_IT", "1000"))
N_CORES <- as.integer(Sys.getenv("SLOWFAST_CORES", max(1L, parallel::detectCores() - 1L)))
NCT_GAMMA <- 0.25
SEED <- 20261019L

dat <- sf_load_data()
st <- ci_setup(sf_symptoms)
X_obs <- as.matrix(dat[, sf_symptoms])
z <- dat$srl_z
high <- dat$high_srl
n_lo <- sum(high == 0); n_hi <- sum(high == 1)
out_dir <- sf_out("step2_simulation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------
# Generating conditions
# ------------------------------------------------------------------------

Jvec <- function(M) M[st$pair_index]
comp_null <- sf_read_perspective("equal_J_equal_h")
comp_field <- sf_read_perspective("equal_J_free_h")
comp_free <- sf_read_perspective("free_J_free_h")

Z_const <- matrix(1, length(z), 1)
Z_group <- cbind(1, high)
Z_srl <- cbind(1, z)

conditions <- list(
  null = list(B = matrix(c(comp_null$h_low, Jvec(comp_null$J_low)), ncol = 1), Z = Z_const),
  fixed_field = list(B = cbind(c(comp_field$h_low, Jvec(comp_field$J_low)),
                               c(comp_field$h_high - comp_field$h_low, rep(0, st$q))), Z = Z_group),
  true_coupling = list(B = cbind(c(comp_free$h_low, Jvec(comp_free$J_low)),
                                 c(comp_free$h_high - comp_free$h_low,
                                   Jvec(comp_free$J_high) - Jvec(comp_free$J_low))), Z = Z_group)
)

step1_file <- file.path(sf_out("step1_model_fits"), "fits.rds")
if (file.exists(step1_file)) {
  s1 <- readRDS(step1_file)$unadjusted
  conditions$context <- list(B = s1$ctx_item$B, Z = Z_srl)
  conditions$context_coupling <- list(B = s1$ctx_coupling$B, Z = Z_srl)
  param_source <- "step 1 fits (ctx_item, ctx_coupling)"
} else {
  cal <- sf_calibrate_context(z, st)
  conditions$context <- list(B = cal$B, Z = Z_srl)
  param_source <- sprintf("calibrated to perspective-paper estimates (max error %.1e)", cal$max_error)
}
conditions <- conditions[intersect(c("null", "fixed_field", "context", "true_coupling", "context_coupling"),
                                   names(conditions))]
message("Context-condition parameters: ", param_source)


# ------------------------------------------------------------------------
# The perspective-paper analysis applied to one data set (or to expected statistics)
# ------------------------------------------------------------------------

two_group_from_X <- function(X) {
  Tm <- ci_features(X, st$pair_index)
  fits <- ci_two_group_suffstat_fit(
    list(colSums(Tm[high == 0, , drop = FALSE]), colSums(Tm[high == 1, , drop = FALSE])),
    list(n_lo, n_hi), st)
  ci_two_group_summary(fits, st)
}

network_stats <- function(X) {
  # The statistics the perspective paper's NCT used: IsingFit, gamma = .25, AND rule.
  fit_net <- function(D) {
    capture.output(r <- IsingFit::IsingFit(D, gamma = NCT_GAMMA, AND = TRUE,
                                           plot = FALSE, progressbar = FALSE))
    r$weiadj
  }
  W1 <- fit_net(X[high == 0, , drop = FALSE]); W2 <- fit_net(X[high == 1, , drop = FALSE])
  ut <- upper.tri(W1)
  data.frame(M = max(abs(W2[ut] - W1[ut])),
             strength_low = sum(abs(W1[ut])), strength_high = sum(abs(W2[ut])),
             S_signed = sum(abs(W2[ut])) - sum(abs(W1[ut])),
             S = abs(sum(abs(W2[ut])) - sum(abs(W1[ut]))))
}

population_level <- function(cond) {
  s_lo <- ci_expected_stats(cond$B, cond$Z[high == 0, , drop = FALSE], st)
  s_hi <- ci_expected_stats(cond$B, cond$Z[high == 1, , drop = FALSE], st)
  fits <- ci_two_group_suffstat_fit(list(s_lo, s_hi), list(n_lo, n_hi), st)
  out <- ci_two_group_summary(fits, st)
  # At the population level LR is the non-centrality; finite-sample E[LR] ~ df + lambda
  names(out)[names(out) == "LR_J"] <- "lambda_J"
  names(out)[names(out) == "LR_h"] <- "lambda_h"
  out$expected_LR_J <- out$lambda_J + st$q
  out$expected_LR_h <- out$lambda_h + st$p
  out
}


# ------------------------------------------------------------------------
# Observed data
# ------------------------------------------------------------------------

message("Analyzing observed data ...")
observed <- cbind(two_group_from_X(X_obs), network_stats(X_obs))
sf_write_csv(observed, file.path(out_dir, "observed_summary.csv"))
if (sf_mode() == "helius") {
  sf_write_csv(sf_perspective_observed(), file.path(out_dir, "perspective_published_summary.csv"))
}


# ------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------

pop <- do.call(rbind, lapply(names(conditions), function(cn) {
  message("Population-level fit: ", cn)
  cbind(condition = cn, population_level(conditions[[cn]]))
}))
sf_write_csv(pop, file.path(out_dir, "population_level.csv"))

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

run_condition <- function(cn) {
  cond <- conditions[[cn]]
  message(sprintf("Replications: %s (%d reps, %d with network statistics)",
                  cn, N_REPS, N_NETWORK_REPS))
  res <- parallel::mclapply(seq_len(N_REPS), function(r) {
    X <- ci_sample(cond$B, cond$Z, st)
    out <- two_group_from_X(X)
    if (r <= N_NETWORK_REPS) out <- cbind(out, network_stats(X))
    else out <- cbind(out, data.frame(M = NA, strength_low = NA, strength_high = NA,
                                      S_signed = NA, S = NA))
    cbind(condition = cn, rep = r, out)
  }, mc.cores = N_CORES, mc.set.seed = TRUE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) warning(cn, ": ", sum(bad), " replications failed")
  do.call(rbind, res[!bad])
}

reps <- do.call(rbind, lapply(names(conditions), run_condition))
sf_write_csv(reps, file.path(out_dir, "replications.csv"))


# ------------------------------------------------------------------------
# Compare observed with each condition
# ------------------------------------------------------------------------

numeric_metrics <- c("LR_J", "LR_h", "dBIC_commonJ_vs_free", "dBIC_commonJ_vs_commonh",
                     "sum_dJ_free", "n_dh_free_positive", "n_dh_commonJ_positive",
                     "dh_slp_free", "dh_slp_commonJ", "M", "S_signed")

# Null-condition benchmarks for the network statistics (parametric NCT)
null_reps <- reps[reps$condition == "null", ]
crit_M <- stats::quantile(null_reps$M, .95, na.rm = TRUE)
crit_S <- stats::quantile(abs(null_reps$S_signed), .95, na.rm = TRUE)

summ <- do.call(rbind, lapply(names(conditions), function(cn) {
  r <- reps[reps$condition == cn, ]
  num <- do.call(rbind, lapply(numeric_metrics, function(v) {
    x <- r[[v]]; x <- x[is.finite(x)]
    data.frame(condition = cn, metric = v, observed = observed[[v]],
               mean = mean(x), q025 = stats::quantile(x, .025), q975 = stats::quantile(x, .975),
               p_sim_ge_obs = mean(x >= observed[[v]]),
               two_sided_p = min(1, 2 * min(mean(x >= observed[[v]]), mean(x <= observed[[v]]))),
               stringsAsFactors = FALSE)
  }))
  cat_rows <- data.frame(
    condition = cn,
    metric = c("P(BIC picks observed model)", "P(AIC picks observed model)",
               "P(M exceeds null 95%)", "P(|S| exceeds null 95%)"),
    observed = c(observed$bic_best, observed$aic_best,
                 as.character(observed$M > crit_M), as.character(abs(observed$S_signed) > crit_S)),
    mean = c(mean(r$bic_best == observed$bic_best), mean(r$aic_best == observed$aic_best),
             mean(r$M > crit_M, na.rm = TRUE), mean(abs(r$S_signed) > crit_S, na.rm = TRUE)),
    q025 = NA, q975 = NA, p_sim_ge_obs = NA, two_sided_p = NA, stringsAsFactors = FALSE)
  rbind(num, cat_rows)
}))
rownames(summ) <- NULL
sf_write_csv(summ, file.path(out_dir, "conditions_vs_observed.csv"))


# ------------------------------------------------------------------------
# Optional: full permutation NCT on a few replications
# ------------------------------------------------------------------------

if (N_FULL_NCT_REPS > 0L) {
  nct_rows <- list()
  for (cn in names(conditions)) for (r in seq_len(N_FULL_NCT_REPS)) {
    message(sprintf("Full NCT: %s rep %d", cn, r))
    X <- ci_sample(conditions[[cn]]$B, conditions[[cn]]$Z, st)
    res <- NetworkComparisonTest::NCT(
      data1 = X[high == 0, ], data2 = X[high == 1, ], gamma = NCT_GAMMA, it = FULL_NCT_IT,
      binary.data = TRUE, paired = FALSE, weighted = TRUE, AND = TRUE, abs = TRUE,
      test.edges = FALSE, progressbar = FALSE, make.positive.definite = TRUE)
    nct_rows[[length(nct_rows) + 1L]] <- data.frame(
      condition = cn, rep = r, M = res$nwinv.real, M_p = res$nwinv.pval,
      S = res$glstrinv.real, S_p = res$glstrinv.pval)
  }
  sf_write_csv(do.call(rbind, nct_rows), file.path(out_dir, "full_nct.csv"))
}


# Figure: see R/revision_2026/figures/fig_helius_mechanisms.R


writeLines(c(
  paste("Data mode:", sf_mode()),
  paste("Context parameters:", param_source),
  paste("Replications per condition:", N_REPS),
  paste("Replications with network statistics:", N_NETWORK_REPS),
  paste("Full permutation NCT replications:", N_FULL_NCT_REPS)
), file.path(out_dir, "run_info.txt"))
message("Step 2 complete. Results in ", normalizePath(out_dir))
