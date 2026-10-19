# ============================================================
# R/revision_2026/empirical/02_fit_context_models.R
#
# Step 1: Fit the cross-sectional slow-fast model to the data and compare it
# with "rewiring" accounts on the same exact likelihood.
#
# The slow-fast model's cross-sectional (fast-layer equilibrium) prediction:
# symptoms share one interaction matrix J, and each person's slow context
# shifts their activations,
#
#   h_n = h0 + gamma * SRL_n        (ctx_item; ctx_common ties gamma)
#
# Rival accounts let interactions themselves change with context, either
# between median-split groups (as in standard network comparisons) or
# continuously with SRL.
#
# Models (all on the same complete-case sample):
#
#   group_commonJ     High/Low split moderates h only        perspective equal_J_free_h
#   group_free        High/Low split moderates h and J       perspective free_J_free_h
#   ctx_common        SRL shifts all h by one common gamma   paper's homogeneous field
#   ctx_item          SRL shifts h with gamma_i per symptom  core slow-fast model
#   ctx_item_quad     + SRL^2 on h (nonlinear threshold shift)
#   ctx_item_groupJ   ctx_item + High/Low split moderates J  does the group J
#                                                            difference survive?
#   ctx_coupling      SRL moderates h and J continuously     continuous rewiring
#   ctx_quad_coupling ctx_coupling + SRL^2 on h              robustness
#
# Key tests:
#   group_free vs group_commonJ      reproduces the perspective-paper LR test of J
#   ctx_item_groupJ vs ctx_item      the same J test once thresholds track SRL
#   ctx_coupling vs ctx_item         continuous SRL moderation of coupling
#
# Each model is fitted unadjusted and adjusted (age, gender, ethnicity shift
# the activations; they never moderate J).
#
# Run from this repository's root (slow-fast-computational-model):
#   Rscript R/revision_2026/empirical/02_fit_context_models.R
#   SLOWFAST_DATA_MODE=synthetic Rscript R/revision_2026/empirical/02_fit_context_models.R
# ============================================================

source("R/revision_2026/empirical/00_context_ising_core.R")
source("R/revision_2026/empirical/01_data.R")

RUN_CV <- as.logical(Sys.getenv("SLOWFAST_RUN_CV", "TRUE"))
CV_FOLDS <- 5L
SEED <- 20261019L

dat <- sf_load_data()
st <- ci_setup(sf_symptoms)
X <- as.matrix(dat[, sf_symptoms])
T_obs <- ci_features(X, st$pair_index)
n <- nrow(dat)
message("Analytic sample: n = ", n)

cov_data <- data.frame(
  srl = dat$srl_z,
  srl2 = dat$srl_z^2 - mean(dat$srl_z^2),
  high = dat$high_srl,
  age = dat$age_z,
  female = dat$female,
  nondutch = dat$non_dutch
)

model_specs <- list(
  group_commonJ     = list(high = "h_item"),
  group_free        = list(high = "hJ_item"),
  ctx_common        = list(srl = "h_common"),
  ctx_item          = list(srl = "h_item"),
  ctx_item_quad     = list(srl = "h_item", srl2 = "h_item"),
  ctx_item_groupJ   = list(srl = "h_item", high = "J_item"),
  ctx_coupling      = list(srl = "hJ_item"),
  ctx_quad_coupling = list(srl = "hJ_item", srl2 = "h_item")
)
adjust_effects <- list(age = "h_item", female = "h_item", nondutch = "h_item")

build <- function(effects, rows = seq_len(n)) {
  covs <- names(effects)
  S <- ci_structure(st, covs, effects)
  Z <- cbind(1, as.matrix(cov_data[rows, covs, drop = FALSE]))
  colnames(Z) <- c("(Intercept)", covs)
  list(S = S, Z = Z)
}

# Adjusted models need complete covariates; unadjusted models use everyone.
rows_for <- function(adjusted) {
  if (adjusted) which(stats::complete.cases(cov_data[, names(adjust_effects)])) else seq_len(n)
}

fit_all <- function(adjusted) {
  fits <- list()
  rows <- rows_for(adjusted)
  if (adjusted && length(rows) < n) {
    message("Adjusted models: n = ", length(rows), " (", n - length(rows),
            " dropped for missing covariates)")
  }
  for (m in names(model_specs)) {
    eff <- model_specs[[m]]
    if (adjusted) eff <- c(eff, adjust_effects)
    b <- build(eff, rows)
    message(sprintf("Fitting %-18s (%s, %d parameters)", m,
                    if (adjusted) "adjusted" else "unadjusted", attr(b$S, "npar")))
    t0 <- Sys.time()
    fits[[m]] <- ci_fit(b$S, NULL, b$Z, st, label = m, T_obs = T_obs[rows, , drop = FALSE])
    fits[[m]]$Z_cols <- colnames(b$Z)
    message(sprintf("   logLik %.2f   %.1f s", fits[[m]]$logLik,
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  fits
}

lrt_table <- function(f) {
  rbind(
    cbind(test = "J differs between groups (perspective-paper test)", ci_lrt(f$group_free, f$group_commonJ)),
    cbind(test = "J differs between groups, given SRL-shifted thresholds", ci_lrt(f$ctx_item_groupJ, f$ctx_item)),
    cbind(test = "J moderated continuously by SRL", ci_lrt(f$ctx_coupling, f$ctx_item)),
    cbind(test = "J moderated by SRL, given nonlinear thresholds", ci_lrt(f$ctx_quad_coupling, f$ctx_item_quad)),
    cbind(test = "Nonlinear (SRL^2) threshold shift", ci_lrt(f$ctx_item_quad, f$ctx_item)),
    cbind(test = "Symptom-specific vs common gamma", ci_lrt(f$ctx_item, f$ctx_common))
  )
}

# Effects of a covariate on the J block: estimates, SEs, and their sum
coupling_effects <- function(fit, covariate) {
  col <- match(covariate, colnames(fit$S))
  rows <- which(!st$is_h)
  idx <- fit$S[rows, col]
  est <- fit$B[rows, col]
  se <- sqrt(diag(fit$vcov))[idx]
  w <- numeric(fit$npar); w[idx] <- 1
  sum_se <- sqrt(as.numeric(t(w) %*% fit$vcov %*% w))
  list(
    edges = data.frame(edge = st$feat_names[rows], estimate = est, se = se,
                       z = est / se, p = 2 * stats::pnorm(-abs(est / se)),
                       p_holm = stats::p.adjust(2 * stats::pnorm(-abs(est / se)), "holm")),
    total = data.frame(covariate = covariate, sum_estimate = sum(est), sum_se = sum_se,
                       z = sum(est) / sum_se, p = 2 * stats::pnorm(-abs(sum(est) / sum_se)))
  )
}

threshold_effects <- function(fit, covariate) {
  col <- match(covariate, colnames(fit$S))
  rows <- which(st$is_h)
  idx <- fit$S[rows, col]
  data.frame(symptom = st$symptoms, gamma = fit$B[rows, col],
             se = sqrt(diag(fit$vcov))[idx])
}

out_dir <- sf_out("step1_model_fits")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
all_fits <- list()

for (adjusted in c(FALSE, TRUE)) {
  tag <- if (adjusted) "adjusted" else "unadjusted"
  fits <- fit_all(adjusted)
  all_fits[[tag]] <- fits

  comp <- ci_compare(fits)
  if (RUN_CV) {
    message("Cross-validating (", CV_FOLDS, " folds) ...")
    comp$cv_loglik_per_person <- sapply(names(fits), function(m) {
      eff <- model_specs[[m]]; if (adjusted) eff <- c(eff, adjust_effects)
      rows <- rows_for(adjusted)
      b <- build(eff, rows)
      ci_cv_loglik(b$S, X[rows, , drop = FALSE], b$Z, st, K = CV_FOLDS, seed = SEED)
    })
    comp$cv_delta_vs_best <- comp$cv_loglik_per_person - max(comp$cv_loglik_per_person)
  }
  sf_write_csv(comp, file.path(out_dir, paste0("model_comparison_", tag, ".csv")))
  sf_write_csv(lrt_table(fits), file.path(out_dir, paste0("likelihood_ratio_tests_", tag, ".csv")))

  ce_cont <- coupling_effects(fits$ctx_coupling, "srl")
  ce_grp <- coupling_effects(fits$ctx_item_groupJ, "high")
  ce_free <- coupling_effects(fits$group_free, "high")
  sf_write_csv(ce_cont$edges, file.path(out_dir, paste0("coupling_by_srl_edges_", tag, ".csv")))
  sf_write_csv(ce_grp$edges, file.path(out_dir, paste0("coupling_group_given_ctx_edges_", tag, ".csv")))
  sf_write_csv(rbind(
    cbind(model = "group_free (High - Low)", ce_free$total),
    cbind(model = "ctx_item_groupJ (High - Low, thresholds track SRL)", ce_grp$total),
    cbind(model = "ctx_coupling (per SD SRL)", ce_cont$total)
  ), file.path(out_dir, paste0("total_coupling_change_", tag, ".csv")))
  sf_write_csv(threshold_effects(fits$ctx_item, "srl"),
               file.path(out_dir, paste0("gamma_ctx_item_", tag, ".csv")))
}

# Fitted parameters (aggregate, not individual-level) for steps 2 and 3
saveRDS(lapply(all_fits, function(fl) lapply(fl, function(f)
  f[c("label", "S", "par", "B", "vcov", "logLik", "npar", "n", "AIC", "BIC", "Z_cols")])),
  file.path(out_dir, "fits.rds"))

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Step 1 complete. Results in ", normalizePath(out_dir))
