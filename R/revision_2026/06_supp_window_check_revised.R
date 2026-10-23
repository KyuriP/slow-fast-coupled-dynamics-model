# Exploratory observation-window analysis, excluded from the submitted manuscript.
# Use the revised nine-symptom model rather than the earlier twelve-symptom model.
#
# Each person's context follows a mean-reverting process with no shocks or feedback.
# After burn-in, average context over W steps to obtain Pbar.
# Generate the symptom profile conditional on Pbar and compare network estimates
# that omit or include it. This assumes symptoms respond to the window average.
# There are no fixed person-level contextual offsets in this design.
#
# The low diffusion scale is 0.04, matching Simulations 2 and 3.
# The high scale is 0.65, chosen to make the window-length effect more visible.
# Both conditions are retained, including small effects at the low scale.
# The generating coupling matrix stays fixed across replicates.
# Outputs are saved under res/revision_2026/window_check.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

source("R/revision_2026/utils_uncentered01_model.R")   # simulate_fast_sweep(), fit_edges()
source("R/revision_2026/00_parameters_uncentered01.R")  # tau (+1.3 shift), omega, gamma, symptoms, N

# ------------------------------------------------------------------------
# Slow-field OU parameters -- identical to Simulations 2/3 so tau_P means
# the same thing as elsewhere in the paper.
# ------------------------------------------------------------------------
kappa   <- 0.20
dt      <- 0.02
P_base  <- 0
tau_P_steps <- 1 / (kappa * dt)   # = 250 steps, same interpretation as the old check

sigma_P_levels <- c(low = 0.04, high = 0.65)

burn_in_ou_steps <- round(4 * tau_P_steps)   # let P reach OU stationarity first
window_steps     <- c(1, 25, 100, 250, 500, 1000, 2500, 5000)  # W, in steps

T_burn_symptom <- 200L   # fast-layer burn-in per person, same as Sim 1 / Sim 4

n_person <- 2000L        # per (sigma_P, W) condition
n_reps   <- 5L            # Monte Carlo repeats of the same DGP, for SEs

iu <- which(upper.tri(matrix(0, N, N)))
true_omega_edges <- omega[iu]
true_gs <- sum(abs(true_omega_edges))
# All nonzero true couplings are positive (same design as Simulation 4),
# so the signed sum equals the absolute sum here -- true_total_coupling
# is kept as its own name because the summary/plotting code below reports
# signed total coupling as the primary recovery metric, matching
# Simulation 4's fix, while true_gs is kept only as a legacy diagnostic.
true_total_coupling <- sum(true_omega_edges)

n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))
RNGkind("L'Ecuyer-CMRG")
set.seed(20260826)

# ------------------------------------------------------------------------
# Vectorized OU window-average: n_person independent people at once.
# Returns Pbar_i for each person (the observed, window-averaged covariate).
# ------------------------------------------------------------------------
simulate_ou_window_means <- function(sigma_P, n_ind, W_steps, burn_in) {
  P <- rep(P_base, n_ind)
  total_steps <- burn_in + W_steps
  running_sum <- numeric(n_ind)
  for (t in seq_len(total_steps)) {
    P <- P + kappa * (P_base - P) * dt + sigma_P * sqrt(dt) * rnorm(n_ind)
    if (t > burn_in) running_sum <- running_sum + P
  }
  running_sum / W_steps
}

# ------------------------------------------------------------------------
# One person's symptom snapshot: same logic as Simulation 4's
# simulate_person(), just conditioning on Pbar_i instead of a fixed/drawn P.
# ------------------------------------------------------------------------
simulate_person_given_Pbar <- function(Pbar) {
  S <- rbinom(N, size = 1, prob = 0.5)
  for (t in seq_len(T_burn_symptom)) {
    S <- simulate_fast_sweep(S, tau, omega, gamma, Pbar)
  }
  S
}

# ------------------------------------------------------------------------
# One (sigma_P, W, rep) condition: simulate Pbar for n_person people,
# generate their symptom profiles, fit symptom-only vs context-adjusted
# networks, and summarize.
# ------------------------------------------------------------------------
run_condition <- function(sigma_P_name, W_steps, rep) {
  sigma_P <- sigma_P_levels[[sigma_P_name]]
  seed <- 1e6 * which(names(sigma_P_levels) == sigma_P_name) + 1000 * rep + W_steps
  set.seed(seed)

  Pbar <- simulate_ou_window_means(sigma_P, n_person, W_steps, burn_in_ou_steps)

  rows <- parallel::mclapply(seq_len(n_person), function(i) simulate_person_given_Pbar(Pbar[i]),
                              mc.cores = n_cores, mc.set.seed = TRUE)
  failed <- vapply(rows, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed)) stop(sprintf("%d of %d person-simulations failed.", sum(failed), n_person))
  S <- do.call(rbind, rows)

  omega_hat_naive    <- fit_edges(S, Pvec = NULL)
  omega_hat_adjusted <- fit_edges(S, Pvec = Pbar)

  tibble(
    sigma_P_level = sigma_P_name, sigma_P = sigma_P, W_steps = W_steps, rep = rep,
    sd_pbar = sd(Pbar),
    total_naive = sum(omega_hat_naive[iu]),
    total_adjusted = sum(omega_hat_adjusted[iu]),
    spurious_naive = sum(abs(omega_hat_naive[iu][true_omega_edges == 0])),
    spurious_adjusted = sum(abs(omega_hat_adjusted[iu][true_omega_edges == 0])),
    # kept as diagnostics
    gs_naive = sum(abs(omega_hat_naive[iu])),
    gs_adjusted = sum(abs(omega_hat_adjusted[iu])),
    mae_naive_all = mean(abs(omega_hat_naive[iu] - true_omega_edges)),
    mae_adjusted_all = mean(abs(omega_hat_adjusted[iu] - true_omega_edges)),
    mae_naive_true_zero = mean(abs(omega_hat_naive[iu][true_omega_edges == 0])),
    mae_adjusted_true_zero = mean(abs(omega_hat_adjusted[iu][true_omega_edges == 0])),
    n_phantom_naive_gt_01 = sum(true_omega_edges == 0 & abs(omega_hat_naive[iu]) > 0.1),
    n_phantom_adjusted_gt_01 = sum(true_omega_edges == 0 & abs(omega_hat_adjusted[iu]) > 0.1)
  )
}

# ------------------------------------------------------------------------
# Run: 2 sigma_P levels x 8 window lengths x n_reps, sequential over
# conditions (each condition is already parallelized internally over
# n_person via mclapply above).
# ------------------------------------------------------------------------
design <- expand_grid(sigma_P_level = names(sigma_P_levels), W_steps = window_steps, rep = seq_len(n_reps))
cat(sprintf("Running %d conditions (%d sigma_P levels x %d window lengths x %d reps), n_person=%d each...\n",
            nrow(design), length(sigma_P_levels), length(window_steps), n_reps, n_person))

results <- design |>
  purrr::pmap_dfr(function(sigma_P_level, W_steps, rep) {
    cat(sprintf("  sigma_P=%s, W=%d, rep=%d\n", sigma_P_level, W_steps, rep))
    run_condition(sigma_P_level, W_steps, rep)
  })

dir.create("res/revision_2026/window_check", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

saveRDS(list(results = results,
             params = list(kappa = kappa, dt = dt, P_base = P_base,
                            sigma_P_levels = sigma_P_levels,
                            burn_in_ou_steps = burn_in_ou_steps,
                            window_steps = window_steps,
                            T_burn_symptom = T_burn_symptom,
                            n_person = n_person, n_reps = n_reps,
                            true_omega = omega, true_gs = true_gs)),
        "res/revision_2026/window_check/window_check_raw.rds")

summary_tbl <- results |>
  group_by(sigma_P_level, sigma_P, W_steps) |>
  summarise(across(
    c(sd_pbar, total_naive, total_adjusted, spurious_naive, spurious_adjusted,
      gs_naive, gs_adjusted, mae_naive_all, mae_adjusted_all,
      mae_naive_true_zero, mae_adjusted_true_zero,
      n_phantom_naive_gt_01, n_phantom_adjusted_gt_01),
    list(mean = mean, se = ~ sd(.x) / sqrt(length(.x))),
    .names = "{.col}_{.fn}"
  ), .groups = "drop") |>
  mutate(W_over_tauP = W_steps / tau_P_steps,
         gap_total = total_naive_mean - total_adjusted_mean,
         gap_spurious = spurious_naive_mean - spurious_adjusted_mean,
         # kept as diagnostics
         gap_gs = gs_naive_mean - gs_adjusted_mean,
         gap_mae_true_zero = mae_naive_true_zero_mean - mae_adjusted_true_zero_mean)

write.csv(summary_tbl, "res/revision_2026/window_check/window_check_summary.csv", row.names = FALSE)

cat("\n=== SUMMARY (naive vs adjusted, by sigma_P and window length) ===\n")
print(summary_tbl |> select(sigma_P_level, W_steps, W_over_tauP, sd_pbar_mean,
                             total_naive_mean, total_adjusted_mean, gap_total,
                             spurious_naive_mean, spurious_adjusted_mean, gap_spurious))

cat(sprintf("\ntrue total coupling = %.2f\n", true_total_coupling))
cat("\nCheck sd_pbar_mean first: if it doesn't decrease appreciably across\n")
cat("W_steps within a sigma_P level, the window-averaging manipulation isn't\n")
cat("doing anything at that diffusion scale, and gap_total / gap_spurious\n")
cat("shouldn't be expected to shrink either -- a genuine, reportable null\n")
cat("result at sigma_P=0.04 (the scale used in Sims 2-3), not a bug.\n")
cat("Compare against sigma_P=0.65, which should show a visible sd_pbar\n")
cat("decrease if the OU/window mechanics are working.\n")
cat("\nReport whatever gap_total / gap_spurious actually come out as here --\n")
cat("these replace the old absolute-global-strength gap_gs / gap_mae_true_zero\n")
cat("as the primary recovery metrics, for the same reason as Simulation 4:\n")
cat("summing |estimate| over true-zero edges makes symmetric sampling noise\n")
cat("look like one-directional inflation. Do NOT reuse the pre-revision\n")
cat("N=12 script's numbers (gap 12.3->10.8->7.7->3.9, mae 0.132->0.024) --\n")
cat("those belong to a different model entirely.\n")

# ------------------------------------------------------------------------
# Figure: 3-panel summary, one column per sigma_P level, log10(W/tau_P) x-axis
# ------------------------------------------------------------------------
col_raw <- "grey15"
col_adj <- "#1565C0"

theme_pub_local <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank()
    )
}

x_lab <- expression(W / tau[P] ~ "(observation window / context timescale)")

panelA <- summary_tbl |>
  ggplot(aes(W_over_tauP, sd_pbar_mean, colour = sigma_P_level)) +
  geom_line() + geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c(low = "#2E7D32", high = "#C62828"),
                       labels = c(low = "sigma_P = 0.04 (main-text scale)", high = "sigma_P = 0.65")) +
  labs(x = x_lab, y = expression(SD(bar(P)[W])), title = "(A) Between-person SD of window-averaged context", colour = NULL) +
  theme_pub_local()

panelB <- summary_tbl |>
  select(sigma_P_level, W_over_tauP, total_naive_mean, total_adjusted_mean) |>
  pivot_longer(c(total_naive_mean, total_adjusted_mean), names_to = "estimator", values_to = "total_coupling") |>
  mutate(estimator = if_else(estimator == "total_naive_mean", "Symptom-only", "Context-adjusted")) |>
  ggplot(aes(W_over_tauP, total_coupling, colour = estimator, linetype = sigma_P_level)) +
  geom_hline(yintercept = true_total_coupling, colour = "grey50", linetype = "dotted") +
  geom_line() + geom_point(size = 1.8) +
  scale_x_log10() +
  scale_colour_manual(values = c("Symptom-only" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = x_lab, y = "Estimated total coupling", title = "(B) Estimated total coupling", colour = NULL, linetype = "sigma_P") +
  theme_pub_local()

panelC <- summary_tbl |>
  select(sigma_P_level, W_over_tauP, spurious_naive_mean, spurious_adjusted_mean) |>
  pivot_longer(c(spurious_naive_mean, spurious_adjusted_mean), names_to = "estimator", values_to = "spurious") |>
  mutate(estimator = if_else(estimator == "spurious_naive_mean", "Symptom-only", "Context-adjusted")) |>
  ggplot(aes(W_over_tauP, spurious, colour = estimator, linetype = sigma_P_level)) +
  geom_line() + geom_point(size = 1.8) +
  scale_x_log10() +
  scale_colour_manual(values = c("Symptom-only" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = x_lab, y = "Spurious absolute coupling on true-zero edges", title = "(C) Spurious coupling on uncoupled pairs", colour = NULL, linetype = "sigma_P") +
  theme_pub_local()

library(patchwork)
fig_window <- (panelA | panelB | panelC) & theme(legend.position = "bottom")

ggsave("figs/revision_2026/Figure_S_window_check_revised.pdf", fig_window, width = 13.5, height = 4.8)
ggsave("figs/revision_2026/Figure_S_window_check_revised.png", fig_window, width = 13.5, height = 4.8, dpi = 200)

cat("\nDone. Files:\n")
cat("  res/revision_2026/window_check/window_check_raw.rds\n")
cat("  res/revision_2026/window_check/window_check_summary.csv\n")
cat("  figs/revision_2026/Figure_S_window_check_revised.pdf (+ .png)\n")
