# Simulation 3. Compare recovery with feedback off and with b = 0.5.
# The symptom model, perturbation and contextual recovery parameters match Simulation 2.
# The exploratory grid 0.3, 0.5 and 0.7 informed the choice of b = 0.5.
# This value gave a visible recovery delay while activation still declined.
#
# Feedback uses the signed difference b * (m_smooth - m_star), replacing
# an earlier positive-part formulation on 2026-08-27.
# m_smooth is an exponentially weighted mean of the active-symptom fraction.
# m_star is the Simulation 1 middle-context mean, 0.2766389, and also supplies
# its initial value. These choices are illustrative rather than empirically fitted.
# Outputs are saved under res/revision_2026/sim3.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

source("R/revision_2026/utils_uncentered01_model.R")
source("R/revision_2026/00_parameters_uncentered01.R")  # tau (+1.3 shift), omega, gamma, symptoms, N

# ------------------------------------------------------------------------
# Slow-state parameters (kappa, sigma_P, dt, shock: identical to Sim 2)
# ------------------------------------------------------------------------
kappa   <- 0.20
sigma_P <- 0.04
dt      <- 0.02
P_base  <- 0

burn_in_steps    <- 200L
post_shock_steps <- 750L
shock_time       <- burn_in_steps + 1L
shock_magnitude  <- 1.0

total_steps <- burn_in_steps + post_shock_steps
n_chains    <- 1000L   # raised from 200: at 200 chains the b=0.3 peak (first
                        # 20 post-shock steps) came out lower than the b=0
                        # peak (4.41 vs 4.73). At 1000 that gap disappeared
                        # (Monte Carlo noise, not a real effect), so 1000 is
                        # the locked chain count for the main-text comparison.

# ------------------------------------------------------------------------
# Feedback parameters used in the manuscript
# ------------------------------------------------------------------------
m_star       <- 0.2766389   # Sim 1 middle-condition (P=0) equilibrium mean_m
alpha_smooth <- 0.05        # EMA rate; half-life = ln(2)/alpha_smooth steps
                             # ~13.9 steps, ~0.28 time units at dt=0.02 --
                             # fast enough to track the post-shock rise/decay,
                             # slow enough to average out single-sweep noise.

b_values <- c(off = 0, on = 0.5)   # Main-text comparison

# ------------------------------------------------------------------------
# One chain: P_t / m_t trajectory with feedback strength b
# ------------------------------------------------------------------------
run_chain <- function(b) {
  P <- P_base
  S <- rbinom(N, size = 1, prob = 0.5)
  m_smooth <- m_star   # start at the reference level, see header note

  P_trace  <- numeric(total_steps)
  M_trace  <- numeric(total_steps)
  ms_trace <- numeric(total_steps)

  for (t in seq_len(total_steps)) {

    # Apply the acute perturbation before symptoms respond at this step.
    if (t == shock_time) {
      P <- P + shock_magnitude
    }

    # Symptoms respond to the current slow field.
    S <- simulate_fast_sweep(S, tau, omega, gamma, P)

    # Record the current state.
    P_trace[t]  <- P
    M_trace[t]  <- symptom_burden(S)
    ms_trace[t] <- m_smooth

    # Update the smoothed burden signal used for feedback.
    m_t <- M_trace[t] / N
    m_smooth <- m_smooth + alpha_smooth * (m_t - m_smooth)

    # Slow recovery/diffusion + feedback for the next step.
    # Signed feedback (2026-08-27, replaces the earlier positive-part
    # clamp): elevated burden pushes P up, burden below m_star pulls it down.
    P <- P + kappa * (P_base - P) * dt + sigma_P * sqrt(dt) * rnorm(1) +
      b * (m_smooth - m_star) * dt
  }
  list(P = P_trace, M = M_trace, m_smooth = ms_trace)
}

# Run independent chains in parallel on macOS/Linux.
# L'Ecuyer-CMRG and mc.set.seed provide separate random-number streams.
# The fork-based mclapply implementation does not support multiple workers on Windows.

n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026L)

tasks <- expand.grid(cond = names(b_values), chain = seq_len(n_chains),
                      stringsAsFactors = FALSE)
n_tasks <- nrow(tasks)

cat(sprintf("Simulating %d chains x %d steps x %d feedback conditions (b = %s), %d tasks total...\n",
            n_chains, total_steps, length(b_values), paste(b_values, collapse = ", "), n_tasks))

run_task <- function(i) {
  cond <- tasks$cond[i]
  b <- b_values[[cond]]
  out <- run_chain(b)
  tibble(feedback = cond, b = b, chain = tasks$chain[i],
         step = seq_len(total_steps), P = out$P, M = out$M, m_smooth = out$m_smooth)
}

chain_tbls <- parallel::mclapply(seq_len(n_tasks), run_task,
                                  mc.cores = n_cores, mc.set.seed = TRUE)

# mclapply returns try-error objects for tasks that fail in a worker
# instead of stopping the whole run, so check for that before trusting output.
failed <- vapply(chain_tbls, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) {
  stop(sprintf("%d of %d parallel tasks failed. First error: %s",
               sum(failed), n_tasks, attr(chain_tbls[[which(failed)[1]]], "condition")$message))
}

traj <- bind_rows(chain_tbls) |>
  mutate(feedback = factor(feedback, levels = names(b_values)),
         time_since_shock = step - shock_time, m = M / N)

# ------------------------------------------------------------------------
# Save + summarize
# ------------------------------------------------------------------------
dir.create("res/revision_2026/sim3", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

saveRDS(list(traj = traj,
             params = list(kappa = kappa, sigma_P = sigma_P, dt = dt, P_base = P_base,
                            b_values = b_values, m_star = m_star, alpha_smooth = alpha_smooth,
                            shock_time = shock_time, shock_magnitude = shock_magnitude,
                            burn_in_steps = burn_in_steps, post_shock_steps = post_shock_steps)),
        "res/revision_2026/sim3/sim3_raw.rds")

summary_tbl <- traj |>
  group_by(feedback, b, step, time_since_shock) |>
  summarise(mean_P = mean(P), se_P = sd(P) / sqrt(n()),
            mean_M = mean(M), se_M = sd(M) / sqrt(n()),
            mean_m = mean(m), .groups = "drop")
write.csv(summary_tbl, "res/revision_2026/sim3/sim3_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------
# Compare recovery with and without feedback.
# ------------------------------------------------------------------------
cat("\n=== PRE-SHOCK BASELINE (last 50 burn-in steps), by feedback condition ===\n")
print(summary_tbl |> filter(time_since_shock < 0, time_since_shock >= -50) |>
        group_by(feedback) |>
        summarise(mean_P = mean(mean_P), mean_M = mean(mean_M), mean_m = mean(mean_m)))

cat("\n=== PEAK RESPONSE (first 20 steps after shock), by feedback condition ===\n")
print(summary_tbl |> filter(time_since_shock >= 0, time_since_shock < 20) |>
        group_by(feedback) |>
        summarise(peak_P = max(mean_P), peak_M = max(mean_M), peak_m = max(mean_m)))

cat("\n=== END-OF-WINDOW RECOVERY (last 50 steps), by feedback condition ===\n")
print(summary_tbl |> filter(time_since_shock >= post_shock_steps - 50) |>
        group_by(feedback) |>
        summarise(mean_P = mean(mean_P), mean_M = mean(mean_M), mean_m = mean(mean_m)))

cat("\nCompare end-of-window activation with the pre-shock baseline.\n")
cat("This analysis examines delayed recovery, not bistability or hysteresis.\n")
cat("Feedback strength is fixed at b = 0.5, selected in\n")
cat("03b_sim_feedback_grid_supplement.R. It is not recalibrated here.\n")

# ------------------------------------------------------------------------
# Figure: P_t / m_t trajectories, b=off vs b=on (0.5) overlaid -- MAIN FIGURE
# ------------------------------------------------------------------------
plot_window <- summary_tbl |> filter(time_since_shock >= -100, time_since_shock <= post_shock_steps)

pP <- ggplot(plot_window, aes(time_since_shock, mean_P, colour = feedback, fill = feedback)) +
  geom_ribbon(aes(ymin = mean_P - se_P, ymax = mean_P + se_P), alpha = 0.2, colour = NA) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c(off = "grey50", on = "firebrick"),
                       labels = c(off = "b = 0", on = sprintf("b = %.2f", b_values[["on"]]))) +
  scale_fill_manual(values = c(off = "grey50", on = "firebrick"), guide = "none") +
  labs(x = NULL, y = expression(P[t]), colour = "Feedback",
       title = "Simulation 3: symptom-to-context feedback slows recovery") +
  theme_classic(base_size = 12)

pM <- ggplot(plot_window, aes(time_since_shock, mean_m, colour = feedback, fill = feedback)) +
  geom_ribbon(aes(ymin = mean_m - se_M / N, ymax = mean_m + se_M / N), alpha = 0.2, colour = NA) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c(off = "grey50", on = "firebrick"),
                       labels = c(off = "b = 0", on = sprintf("b = %.2f", b_values[["on"]]))) +
  scale_fill_manual(values = c(off = "grey50", on = "firebrick"), guide = "none") +
  labs(x = "Steps since shock (dashed line = shock onset)", y = expression(m[t]), colour = "Feedback") +
  theme_classic(base_size = 12)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
p_combined <- patchwork::wrap_plots(pP, pM, ncol = 1) + patchwork::plot_layout(guides = "collect")

ggsave("figs/revision_2026/fig_sim3_feedback_comparison.pdf", p_combined, width = 8, height = 6.5)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim3/sim3_raw.rds\n")
cat("  res/revision_2026/sim3/sim3_summary.csv\n")
cat("  figs/revision_2026/fig_sim3_feedback_comparison.pdf\n")
