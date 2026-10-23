# Simulation 4. Network estimation with fixed or varying context.
# Each person contributes the final symptom profile after 200 fast sweeps.
# The baseline arm fixes P at zero. The other arms use the same profiles
# with P drawn uniformly from -0.6 to 0.6 and either omit or adjust for P.
# The generating network is omega in 00_parameters_uncentered01.R.
# Results are saved in res/revision_2026/sim4.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

source("R/revision_2026/utils_uncentered01_model.R")   # simulate_fast_sweep(), fit_edges()
source("R/revision_2026/00_parameters_uncentered01.R")  # tau (+1.3 shift), omega, gamma, symptoms, N

# ------------------------------------------------------------------------
# Simulation settings
# ------------------------------------------------------------------------
# The initial n = 2,000 run had substantial estimation error even at fixed P.
# This run uses n = 10,000 with the same unregularized estimator.
T_burn    <- 200L
n_person  <- 10000L

P_fixed    <- 0          # baseline arm: everyone simulated at this P
P_range    <- c(-0.6, 0.6)  # naive/adjusted arms: P_i ~ Uniform(P_range),
                             # matching Sim 1's low/high context range

iu <- which(upper.tri(matrix(0, N, N)))   # upper-triangle indices, reused
                                            # for every N x N matrix below
true_omega_edges <- omega[iu]

# ------------------------------------------------------------------------
# One person: simulate forward T_burn sweeps at a given P, return final
# state only (one cross-sectional observation).
# ------------------------------------------------------------------------
simulate_person <- function(P) {
  S <- rbinom(N, size = 1, prob = 0.5)
  for (t in seq_len(T_burn)) {
    S <- simulate_fast_sweep(S, tau, omega, gamma, P)
  }
  S
}

# ------------------------------------------------------------------------
# Generate independent profiles and retain each person's true P.
# The unadjusted estimator does not receive P.
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))
RNGkind("L'Ecuyer-CMRG")
set.seed(2026L)

generate_arm_data <- function(Pvec) {
  n <- length(Pvec)
  rows <- parallel::mclapply(seq_len(n), function(i) simulate_person(Pvec[i]),
                              mc.cores = n_cores, mc.set.seed = TRUE)
  failed <- vapply(rows, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed)) stop(sprintf("%d of %d person-simulations failed.", sum(failed), n))
  S <- do.call(rbind, rows)
  colnames(S) <- symptoms
  list(S = S, Pvec = Pvec)
}

cat("Simulating baseline arm (single context, P=0)...\n")
data_baseline <- generate_arm_data(rep(P_fixed, n_person))

cat("Simulating naive/adjusted arms (pooled, P ~ Uniform(-0.6, 0.6))...\n")
P_pooled <- runif(n_person, P_range[1], P_range[2])
data_pooled <- generate_arm_data(P_pooled)   # same data feeds both naive and adjusted

# ------------------------------------------------------------------------
# Estimate networks
# ------------------------------------------------------------------------
cat("Fitting networks (nodewise logistic regression, symmetrized)...\n")

omega_hat_baseline <- fit_edges(data_baseline$S, Pvec = NULL)
omega_hat_naive     <- fit_edges(data_pooled$S,   Pvec = NULL)                 # P NOT given to the estimator
omega_hat_adjusted  <- fit_edges(data_pooled$S,   Pvec = data_pooled$Pvec)     # P given to the estimator

# ------------------------------------------------------------------------
# Save + summarize
# ------------------------------------------------------------------------
dir.create("res/revision_2026/sim4", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

pair_names <- outer(symptoms, symptoms, paste, sep = "-")[iu]
sym_i <- symptoms[row(omega)[iu]]
sym_j <- symptoms[col(omega)[iu]]

edge_tbl <- bind_rows(
  tibble(arm = "baseline", symptom_i = sym_i, symptom_j = sym_j,
         true_omega = true_omega_edges, est_omega = omega_hat_baseline[iu]),
  tibble(arm = "naive", symptom_i = sym_i, symptom_j = sym_j,
         true_omega = true_omega_edges, est_omega = omega_hat_naive[iu]),
  tibble(arm = "adjusted", symptom_i = sym_i, symptom_j = sym_j,
         true_omega = true_omega_edges, est_omega = omega_hat_adjusted[iu])
) |>
  mutate(arm = factor(arm, levels = c("baseline", "naive", "adjusted")),
         true_nonzero = true_omega != 0,
         error = est_omega - true_omega)

write.csv(edge_tbl, "res/revision_2026/sim4/sim4_edge_estimates.csv", row.names = FALSE)

saveRDS(list(edge_tbl = edge_tbl,
             omega_hat = list(baseline = omega_hat_baseline, naive = omega_hat_naive,
                               adjusted = omega_hat_adjusted),
             true_omega = omega,
             data = list(baseline = data_baseline, pooled = data_pooled),
             params = list(T_burn = T_burn, n_person = n_person,
                            P_fixed = P_fixed, P_range = P_range)),
        "res/revision_2026/sim4/sim4_raw.rds")

summary_tbl <- edge_tbl |>
  group_by(arm) |>
  summarise(
    total_coupling_true = sum(true_omega),
    total_coupling_est  = sum(est_omega),

    spurious_abs_coupling =
      sum(abs(est_omega[!true_nonzero])),

    mae_all = mean(abs(error)),
    mae_true_nonzero = mean(abs(error[true_nonzero])),
    mae_true_zero = mean(abs(error[!true_nonzero])),

    # retained only as diagnostics
    global_strength_est = sum(abs(est_omega)),
    n_phantom_edges_gt_01 =
      sum(!true_nonzero & abs(est_omega) > 0.1),

    .groups = "drop"
  )
write.csv(summary_tbl, "res/revision_2026/sim4/sim4_summary.csv", row.names = FALSE)

cat("\n=== SIMULATION 4 SUMMARY (by arm) ===\n")
print(summary_tbl)

cat("\nCompare 'naive' against 'baseline' and 'adjusted': if pooling across\n")
cat("context without conditioning on P inflates global_strength_est and/or\n")
cat("mae_true_zero (spurious nonzero edges among truly-zero pairs) relative\n")
cat("to both baseline and adjusted, that's the confounding effect this\n")
cat("simulation is designed to show. If 'naive' looks similar to 'baseline'/\n")
cat("'adjusted', the P_range spread may be too small to generate visible\n")
cat("confounding -- widen P_range and rerun. If 'adjusted' doesn't track\n")
cat("'baseline' closely, something is wrong with the P-conditioning itself\n")
cat("(not just a calibration issue).\n")

# ------------------------------------------------------------------------
# Figure: estimated vs true edge weight, one panel per arm
# ------------------------------------------------------------------------
p1 <- ggplot(edge_tbl, aes(true_omega, est_omega, colour = true_nonzero)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                       labels = c(`FALSE` = "true edge = 0", `TRUE` = "true edge != 0"),
                       name = NULL) +
  facet_wrap(~arm, nrow = 1) +
  labs(x = expression(paste("True ", omega[ij])), y = expression(paste("Estimated ", hat(omega)[ij])),
       title = "Simulation 4: omitted-context confounding in estimated symptom networks",
       subtitle = "Dashed line = perfect recovery. Points off the line among grey (true-zero) pairs are spurious edges.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figs/revision_2026/fig_sim4_network_estimation.pdf", p1, width = 11, height = 4.5)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim4/sim4_raw.rds\n")
cat("  res/revision_2026/sim4/sim4_edge_estimates.csv\n")
cat("  res/revision_2026/sim4/sim4_summary.csv\n")
cat("  figs/revision_2026/fig_sim4_network_estimation.pdf\n")
