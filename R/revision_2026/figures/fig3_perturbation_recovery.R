# Simulation 2 trajectories with feedback off.
# Panel A shows contextual recovery. Panel B shows the distribution of
# symptom counts across chains, with the mean trajectory overlaid.
# Counts use the 0 to 9 scale reported in the manuscript.
# Output is Figure3_perturbation_recovery.pdf and .png.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

source("R/revision_2026/figures/theme_publication.R")  # theme_pub(), panel_title(), col_P, col_M

# ------------------------------------------------------------------------
# 0. Load data (Simulation 2 only -- feedback off, b=0)
# ------------------------------------------------------------------------
sim2 <- read_csv("res/revision_2026/sim2/sim2_summary.csv", show_col_types = FALSE)

plot_window <- c(-100, 750)
sim2_w <- sim2 |> filter(time_since_shock >= plot_window[1], time_since_shock <= plot_window[2])

base_sz <- 11

# Smooth the symptom-count mean for display with a two-sided moving average.
# Edge values without a full window are omitted.
roll_mean <- function(x, k = 15) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))

sim2_w <- sim2_w |>
  arrange(time_since_shock) |>
  mutate(mean_M_smooth = roll_mean(mean_M), se_M_smooth = roll_mean(se_M))

# Select a reproducible subset of chains for the contextual trajectories.
set.seed(3)
n_sample_chains <- 25L

sim2_raw <- readRDS("res/revision_2026/sim2/sim2_raw.rds")$traj |>
  filter(time_since_shock >= plot_window[1], time_since_shock <= plot_window[2])
sim2_chains_w <- sim2_raw |> filter(chain %in% sample(unique(sim2_raw$chain), n_sample_chains))

# Bin symptom counts over the same 15-step window used for smoothing.
# Show the distribution across chains and times within each bin.
time_bin_width <- 15

heat_df <- sim2_raw |>
  mutate(time_bin = round(time_since_shock / time_bin_width) * time_bin_width) |>
  count(time_bin, M, name = "n_chains") |>
  group_by(time_bin) |>
  mutate(prop = n_chains / sum(n_chains)) |>
  ungroup()

# ------------------------------------------------------------------------
# 1. y-limits and interval settings
# ------------------------------------------------------------------------
interval_mult <- 1.96  # 95% intervals

range_pad <- function(x, pad = 0.06, lower = -Inf, upper = Inf) {
  r <- range(x, na.rm = TRUE)
  d <- diff(r)
  if (d == 0) d <- 1
  c(max(lower, r[1] - pad * d), min(upper, r[2] + pad * d))
}

p_ylim <- range_pad(
  c(sim2_w$mean_P - interval_mult * sim2_w$se_P, sim2_w$mean_P + interval_mult * sim2_w$se_P),
  pad = 0.06
)

p_ref <- mean(sim2_w$mean_P[sim2_w$time_since_shock < 0], na.rm = TRUE)
# M_ylim/M_ref (ribbon range + pre-perturbation reference line) aren't
# needed -- panel B uses the heatmap below instead, with a fixed 0-9 axis.

# ------------------------------------------------------------------------
# 2. Panel A: slow field P_t -- bold mean + light ribbon + a small
#    faint sample of individual chains
# ------------------------------------------------------------------------
pA <- ggplot(sim2_w, aes(x = time_since_shock, y = mean_P)) +
  geom_hline(yintercept = p_ref, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = sim2_chains_w,
    aes(x = time_since_shock, y = P, group = chain),
    inherit.aes = FALSE,
    colour = col_P, alpha = 0.22, linewidth = 0.18
  ) +
  geom_ribbon(
    aes(ymin = mean_P - interval_mult * se_P, ymax = mean_P + interval_mult * se_P),
    fill = col_P, alpha = 0.12
  ) +
  geom_line(colour = col_P, linewidth = 0.55) +
  coord_cartesian(ylim = p_ylim, clip = "off") +
  labs(title = panel_title("A", "Slow context rises and recovers"), x = NULL, y = "Slow field (P)") +
  theme_pub(base_size = base_sz) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel B. Draw the symptom-count heatmap behind the mean trajectory.
pB <- ggplot(sim2_w, aes(x = time_since_shock, y = mean_M_smooth)) +
  geom_tile(
    data = heat_df,
    aes(x = time_bin, y = M, fill = prop),
    inherit.aes = FALSE,
    width = time_bin_width, height = 1
  ) +
  scale_fill_gradient(
    low = "white", high = col_M, name = "Proportion of chains",
    labels = scales::percent_format(accuracy = 1),
    breaks = scales::breaks_width(0.1),
    guide = guide_colorbar(
      barwidth = unit(5.5, "cm"), barheight = unit(0.42, "cm"),
      title.position = "top", title.hjust = 0.5
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.45) +
  geom_line(colour = "white", linewidth = 1.5, lineend = "round") +
  geom_line(colour = "grey15", linewidth = 0.8) +
  scale_y_continuous(breaks = 0:9, limits = c(-0.5, 9.5), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = panel_title("B", "Symptom activation shifts upward and recovers"),
       x = "Steps since perturbation", y = "Number of active symptoms") +
  theme_pub(base_size = base_sz)

# Stack the two panels with a shared perturbation-relative time axis.
fig_perturbation_recovery <- pA / pB + plot_layout(heights = c(0.9, 1.1))

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave("figs/revision_2026/Figure3_perturbation_recovery.pdf", fig_perturbation_recovery, width = 6.4, height = 6.4)
ggsave("figs/revision_2026/Figure3_perturbation_recovery.png", fig_perturbation_recovery, width = 6.4, height = 6.4, dpi = 300)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure3_perturbation_recovery.pdf (+ .png)\n")
cat("\nThis is a NEW figure, inserted before the existing Figure 3 (regimes,\n")
cat("fig3_regimes.R / Figure3_regimes.pdf) in document order -- that figure\n")
cat("and Figure 4 (network trio) renumber automatically via LaTeX's \\ref,\n")
cat("no manuscript label changes needed beyond adding this figure\n")
cat("environment. See chat for the LaTeX block and the Results paragraph to\n")
cat("insert in the Simulation 2 section.\n")
