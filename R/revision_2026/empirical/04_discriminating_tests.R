# ============================================================
# R/revision_2026/empirical/04_discriminating_tests.R
#
# Step 3: Tests that separate the slow-fast (context) account from a group
# "rewiring" account. Both can fit the perspective paper's single median split;
# they make different predictions elsewhere.
#
# (a) Different splits of SRL. If the apparent J difference arises from
#     unmodelled context variation within groups, its size should change in a
#     predictable way as the split moves (median, terciles, quartiles,
#     quintiles: bottom vs top group). The rewiring account fitted to the
#     median split predicts that people keep the J of their median-split
#     group, whatever the split.
#
# (b) Out-of-sample groupings. The context model is fitted with SRL only. It
#     then predicts, with no new parameters, what the perspective-paper analysis
#     should find for Younger vs Older and Dutch vs Non-Dutch, given how SRL
#     is distributed in those groups.
#
# (c) Pattern of edge differences. For the median split, the context model
#     predicts a specific vector of 36 apparent J differences. We test the
#     observed differences against that prediction (Wald chi-square with the
#     observed covariance), alongside the test against zero.
#
# Predictions come from expected sufficient statistics (exact, no noise);
# prediction intervals for the observed statistic come from replications
# (N_REPS) drawn from each fitted model.
#
# Requires step 1 (02_fit_context_models.R) for the same data mode.
# ============================================================

source("R/revision_2026/empirical/00_context_ising_core.R")
source("R/revision_2026/empirical/01_data.R")

N_REPS <- as.integer(Sys.getenv("SLOWFAST_N_REPS_STEP3", "100"))
N_CORES <- as.integer(Sys.getenv("SLOWFAST_CORES", max(1L, parallel::detectCores() - 1L)))
SEED <- 20261019L

dat <- sf_load_data()
st <- ci_setup(sf_symptoms)
X_obs <- as.matrix(dat[, sf_symptoms])
T_obs <- ci_features(X_obs, st$pair_index)
z <- dat$srl_z
out_dir <- sf_out("step3_discriminating_tests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step1_file <- file.path(sf_out("step1_model_fits"), "fits.rds")
if (!file.exists(step1_file)) stop("Run R/revision_2026/empirical/02_fit_context_models.R first (same data mode).")
s1 <- readRDS(step1_file)$unadjusted

# Generating models and the design each needs for an arbitrary set of people
models <- list(
  context = list(B = s1$ctx_item$B, design = function(idx) cbind(1, z[idx])),
  context_coupling = list(B = s1$ctx_coupling$B, design = function(idx) cbind(1, z[idx])),
  group_rewiring = list(B = s1$group_free$B, design = function(idx) cbind(1, dat$high_srl[idx]))
)


# ------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------

analyse_groups <- function(T_mat, g1, g2) {
  fits <- ci_two_group_suffstat_fit(
    list(colSums(T_mat[g1, , drop = FALSE]), colSums(T_mat[g2, , drop = FALSE])),
    list(length(g1), length(g2)), st)
  list(summary = ci_two_group_summary(fits, st), fits = fits)
}

predict_groups <- function(model, g1, g2) {
  s1e <- ci_expected_stats(model$B, model$design(g1), st)
  s2e <- ci_expected_stats(model$B, model$design(g2), st)
  fits <- ci_two_group_suffstat_fit(list(s1e, s2e), list(length(g1), length(g2)), st)
  out <- ci_two_group_summary(fits, st)
  out$expected_LR_J <- out$LR_J + st$q
  list(summary = out, fits = fits)
}

simulate_interval <- function(model, g1, g2) {
  idx <- c(g1, g2)
  Zs <- model$design(idx)
  res <- parallel::mclapply(seq_len(N_REPS), function(r) {
    X <- ci_sample(model$B, Zs, st)
    Tm <- ci_features(X, st$pair_index)
    analyse_groups(Tm, seq_along(g1), length(g1) + seq_along(g2))$summary[, c("LR_J", "sum_dJ_free")]
  }, mc.cores = N_CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res[!vapply(res, inherits, logical(1), "try-error")])
  data.frame(LR_J_q025 = stats::quantile(res$LR_J, .025), LR_J_q975 = stats::quantile(res$LR_J, .975),
             sum_dJ_q025 = stats::quantile(res$sum_dJ_free, .025),
             sum_dJ_q975 = stats::quantile(res$sum_dJ_free, .975))
}

compare_split <- function(label, g1, g2) {
  message("  ", label, sprintf(" (n = %d vs %d)", length(g1), length(g2)))
  obs <- analyse_groups(T_obs, g1, g2)$summary
  rows <- list(cbind(split = label, source = "observed", obs[, c("LR_J", "sum_dJ_free", "bic_best")],
                     expected_LR_J = NA, LR_J_q025 = NA, LR_J_q975 = NA,
                     sum_dJ_q025 = NA, sum_dJ_q975 = NA,
                     var_srl_g1 = stats::var(z[g1]), var_srl_g2 = stats::var(z[g2])))
  for (mn in names(models)) {
    pr <- predict_groups(models[[mn]], g1, g2)$summary
    iv <- simulate_interval(models[[mn]], g1, g2)
    rows[[length(rows) + 1L]] <- cbind(
      split = label, source = mn,
      data.frame(LR_J = pr$LR_J, sum_dJ_free = pr$sum_dJ_free, bic_best = pr$bic_best),
      expected_LR_J = pr$expected_LR_J, iv,
      var_srl_g1 = stats::var(z[g1]), var_srl_g2 = stats::var(z[g2]))
  }
  out <- do.call(rbind, rows)
  names(out)[names(out) == "LR_J"] <- "LR_J_or_noncentrality"
  out
}

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)


# ------------------------------------------------------------------------
# (a) Different splits of SRL
# ------------------------------------------------------------------------

message("(a) SRL splits")
qz <- function(p) stats::quantile(z, p)
split_defs <- list(
  "median (bottom 1/2 vs top 1/2)" = list(which(dat$high_srl == 0), which(dat$high_srl == 1)),
  "terciles (bottom 1/3 vs top 1/3)" = list(which(z <= qz(1/3)), which(z > qz(2/3))),
  "quartiles (bottom 1/4 vs top 1/4)" = list(which(z <= qz(1/4)), which(z > qz(3/4))),
  "quintiles (bottom 1/5 vs top 1/5)" = list(which(z <= qz(1/5)), which(z > qz(4/5)))
)
splits <- do.call(rbind, lapply(names(split_defs), function(s)
  compare_split(s, split_defs[[s]][[1]], split_defs[[s]][[2]])))
sf_write_csv(splits, file.path(out_dir, "a_srl_splits.csv"))
message("  Note: the median and quantile splits use different n, so compare each ",
        "observed value with its own predictions and intervals, not across rows.")


# ------------------------------------------------------------------------
# (b) Out-of-sample groupings
# ------------------------------------------------------------------------

message("(b) Out-of-sample groupings")
group_defs <- list(
  "Younger vs Older" = list(which(dat$older == 0), which(dat$older == 1)),
  "Dutch vs Non-Dutch" = list(which(dat$non_dutch == 0), which(dat$non_dutch == 1))
)
grouping <- do.call(rbind, lapply(names(group_defs), function(s)
  compare_split(s, group_defs[[s]][[1]], group_defs[[s]][[2]])))
sf_write_csv(grouping, file.path(out_dir, "b_out_of_sample_groupings.csv"))

# Symptom-level activation differences: observed vs predicted from SRL alone
act_rows <- list()
for (s in names(group_defs)) {
  g1 <- group_defs[[s]][[1]]; g2 <- group_defs[[s]][[2]]
  ob <- analyse_groups(T_obs, g1, g2)$fits$equal_J_free_h
  pr <- predict_groups(models$context, g1, g2)$fits$equal_J_free_h
  act_rows[[s]] <- data.frame(grouping = s, symptom = st$symptoms,
                              dh_observed = ob$h2 - ob$h1, dh_predicted_from_srl = pr$h2 - pr$h1)
}
sf_write_csv(do.call(rbind, act_rows), file.path(out_dir, "b_activation_differences.csv"))


# ------------------------------------------------------------------------
# (c) Pattern of edge differences for the median split
# ------------------------------------------------------------------------

message("(c) Edge-difference pattern")
gf <- s1$group_free
J_rows <- which(!st$is_h)
col_high <- match("high", colnames(gf$S))
dJ_obs <- gf$B[J_rows, col_high]
V <- gf$vcov[gf$S[J_rows, col_high], gf$S[J_rows, col_high]]
Vinv <- solve(V)
wald <- function(target) {
  d <- dJ_obs - target
  stat <- as.numeric(t(d) %*% Vinv %*% d)
  c(chisq = stat, df = length(d), p = stats::pchisq(stat, length(d), lower.tail = FALSE))
}
g1 <- which(dat$high_srl == 0); g2 <- which(dat$high_srl == 1)
pred_ctx <- predict_groups(models$context, g1, g2)$fits$free_J_free_h
pred_cc <- predict_groups(models$context_coupling, g1, g2)$fits$free_J_free_h
dJ_ctx <- pred_ctx$J2 - pred_ctx$J1
dJ_cc <- pred_cc$J2 - pred_cc$J1
edge_tab <- data.frame(edge = st$feat_names[J_rows], dJ_observed = dJ_obs,
                       se = sqrt(diag(V)), dJ_pred_context = dJ_ctx,
                       dJ_pred_context_coupling = dJ_cc)
sf_write_csv(edge_tab, file.path(out_dir, "c_edge_differences.csv"))
wtab <- rbind(
  data.frame(target = "zero (no difference)", t(wald(rep(0, st$q)))),
  data.frame(target = "context model prediction", t(wald(dJ_ctx))),
  data.frame(target = "context + coupling model prediction", t(wald(dJ_cc)))
)
wtab$cor_observed_vs_target <- c(NA, stats::cor(dJ_obs, dJ_ctx), stats::cor(dJ_obs, dJ_cc))
sf_write_csv(wtab, file.path(out_dir, "c_edge_pattern_wald.csv"))


# ------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  both <- rbind(cbind(panel = "SRL splits", splits), cbind(panel = "Other groupings", grouping))
  pr <- both[both$source != "observed", ]
  ob <- both[both$source == "observed", ]
  g <- ggplot(pr, aes(y = split, x = (LR_J_q025 + LR_J_q975) / 2, colour = source)) +
    geom_errorbarh(aes(xmin = LR_J_q025, xmax = LR_J_q975), height = .3,
                   position = position_dodge(width = .6)) +
    geom_point(data = ob, aes(x = LR_J_or_noncentrality, y = split), colour = "black",
               shape = 4, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ panel, scales = "free_y", ncol = 1) +
    labs(x = "Likelihood-ratio statistic for equal J (df = 36)", y = NULL,
         colour = "95% interval under",
         caption = "Cross: observed. Bars: replications from each fitted model.") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "fig_discriminating_tests.pdf"), g, width = 8.5, height = 6)
  ggsave(file.path(out_dir, "fig_discriminating_tests.png"), g, width = 8.5, height = 6, dpi = 200)
}

message("Step 3 complete. Results in ", normalizePath(out_dir))
