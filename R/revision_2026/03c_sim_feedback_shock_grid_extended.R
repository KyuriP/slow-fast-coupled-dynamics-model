# Simulation 3. Recovery from a single contextual shock across feedback strengths.
# Uses the same symptom model and shock size as 03_sim_feedback.R.
# Outputs are saved in res/revision_2026/sim3c.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(parallel)
})

source("R/revision_2026/utils_uncentered01_model.R")
source("R/revision_2026/00_parameters_uncentered01.R")  # tau (+1.3 shift), omega, gamma, symptoms, N

# ------------------------------------------------------------------------
# Slow-state parameters -- identical to Sim 2/3
# ------------------------------------------------------------------------
kappa   <- 0.20
sigma_P <- 0.04
dt      <- 0.02
P_base  <- 0

burn_in_steps    <- 200L
# Extended from 1,500 to 2,500 steps on 2026-08-25 because b = 1.25
# had not settled by step 1,500. End-of-window changes are reported below.
post_shock_steps <- 2500L
shock_time       <- burn_in_steps + 1L
shock_magnitude  <- 1.0

total_steps <- burn_in_steps + post_shock_steps
n_chains    <- 300L

m_star       <- 0.2766389   # Sim 1 middle-condition (P=0) equilibrium mean_m
alpha_smooth <- 0.05

# Added b = 0.90, 1.10 and 1.30 on 2026-08-25 to examine the transition.
# Retain the original grid values in the saved results.
b_grid <- c(b000 = 0, b050 = 0.50, b075 = 0.75, b090 = 0.90, b100 = 1.00,
            b110 = 1.10, b125 = 1.25, b130 = 1.30, b150 = 1.50)

# ------------------------------------------------------------------------
# One chain: identical shock-and-recovery design to Sim 2/3
# ------------------------------------------------------------------------
run_chain <- function(b) {
  P <- P_base
  S <- rbinom(N, size = 1, prob = 0.5)
  m_smooth <- m_star

  P_trace <- numeric(total_steps)
  M_trace <- numeric(total_steps)

  for (t in seq_len(total_steps)) {
    if (t == shock_time) {
      P <- P + shock_magnitude
    }

    S <- simulate_fast_sweep(S, tau, omega, gamma, P)

    P_trace[t] <- P
    M_trace[t] <- symptom_burden(S)

    m_t <- M_trace[t] / N
    m_smooth <- m_smooth + alpha_smooth * (m_t - m_smooth)

    P <- P + kappa * (P_base - P) * dt + sigma_P * sqrt(dt) * rnorm(1) +
      b * (m_smooth - m_star) * dt
  }
  list(P = P_trace, M = M_trace)
}

# ------------------------------------------------------------------------
# Run all (b, chain) combinations in parallel -- same mclapply pattern as
# Sim 3 / Sim 4 (fork-based, macOS/Linux only)
# ------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))

RNGkind("L'Ecuyer-CMRG")
set.seed(20260825)

tasks <- expand.grid(cond = names(b_grid), chain = seq_len(n_chains), stringsAsFactors = FALSE)
n_tasks <- nrow(tasks)

cat(sprintf("Simulating %d chains x %d steps x %d feedback levels (b = %s), %d tasks total...\n",
            n_chains, total_steps, length(b_grid), paste(b_grid, collapse = ", "), n_tasks))

run_task <- function(i) {
  cond <- tasks$cond[i]
  b <- b_grid[[cond]]
  out <- run_chain(b)
  tibble(b_name = cond, b = b, chain = tasks$chain[i],
         step = seq_len(total_steps), P = out$P, M = out$M)
}

chain_tbls <- parallel::mclapply(seq_len(n_tasks), run_task, mc.cores = n_cores, mc.set.seed = TRUE)

failed <- vapply(chain_tbls, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) {
  stop(sprintf("%d of %d parallel tasks failed. First error: %s",
               sum(failed), n_tasks, attr(chain_tbls[[which(failed)[1]]], "condition")$message))
}

traj <- bind_rows(chain_tbls) |>
  mutate(b_name = factor(b_name, levels = names(b_grid)),
         time_since_shock = step - shock_time, m = M / N)

dir.create("res/revision_2026/sim3c", recursive = TRUE, showWarnings = FALSE)

saveRDS(list(traj = traj,
             params = list(kappa = kappa, sigma_P = sigma_P, dt = dt, P_base = P_base,
                            b_grid = b_grid, m_star = m_star, alpha_smooth = alpha_smooth,
                            shock_time = shock_time, shock_magnitude = shock_magnitude,
                            burn_in_steps = burn_in_steps, post_shock_steps = post_shock_steps)),
        "res/revision_2026/sim3c/sim3c_raw.rds")

summary_tbl <- traj |>
  group_by(b_name, b, step, time_since_shock) |>
  summarise(mean_P = mean(P), se_P = sd(P) / sqrt(n()),
            mean_M = mean(M), se_M = sd(M) / sqrt(n()),
            mean_m = mean(M / 9), .groups = "drop")
write.csv(summary_tbl, "res/revision_2026/sim3c/sim3c_summary.csv", row.names = FALSE)

cat("\n=== END-OF-WINDOW (last 50 steps), by b ===\n")
print(summary_tbl |> filter(time_since_shock >= post_shock_steps - 50) |>
        group_by(b) |>
        summarise(mean_P = mean(mean_P), mean_m = mean(mean_m)))

cat("\nCompare end-of-window mean_m/mean_P across b: values near the\n")
cat("pre-shock baseline mean recovery; values still clearly elevated after\n")
cat("1500 steps mean a persistently elevated (non-recovering) regime at\n")
cat("that b.\n")

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim3c/sim3c_raw.rds\n")
cat("  res/revision_2026/sim3c/sim3c_summary.csv\n")
