# Feedback, recovery and dependence on the initial state.
# Panel A shows shock recovery across the feedback grid from sim3c.
# Panel B shows the history index from 05_supp_regime_history_dependence.R.
# These use the revised heterogeneous nine-symptom model.
# Output is Figure3_regimes.pdf and .png.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# Use the shared plotting theme.
source("R/revision_2026/figures/theme_publication.R")  # theme_pub(), panel_title(), col_*, pal_feedback

roll_mean <- function(x, k = 15) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
base_sz <- base_size_panel_title  # 10.5, from theme_publication.R's global tokens
interval_mult <- 1.96

# Darker gray keeps the feedback-off trajectory visible under the ribbons.
gray_off_dark <- "#4D4D4D"

# Text sizes for this figure.
axis_text_sz  <- 12
axis_title_sz <- 13
panel_title_sz <- 15
legend_text_sz <- 12.5

# Panel A. Show recovery at b = 0, 0.5, 0.75, 0.9, 1 and 1.3.
# Read the extended 2,500-step runs from sim3c.
# The numerical off/on comparison in the text uses the separate sim3 run,
# so its exact values differ from these trajectories because the draws differ.
# The displayed grid is a subset of the saved grid. The earlier plotting choice
# omitted b = 1.1 because it was still changing near the end of the window.
# A visually flattening trajectory is not proof of asymptotic stability.
display_b_grid <- c(0, 0.50, 0.75, 0.90, 1.00, 1.30)

sim3c <- read_csv("res/revision_2026/sim3c/sim3c_summary.csv", show_col_types = FALSE)
sim3c_raw <- readRDS("res/revision_2026/sim3c/sim3c_raw.rds")$traj

plot_window_grid <- c(-100, 2500)
N <- 9

grid_w <- sim3c |> filter(b %in% display_b_grid, time_since_shock >= plot_window_grid[1], time_since_shock <= plot_window_grid[2]) |>
  arrange(b, time_since_shock) |>
  group_by(b) |>
  mutate(mean_m_smooth = roll_mean(mean_m), se_m_smooth = roll_mean(se_M / N)) |>
  ungroup() |>
  mutate(mean_M_smooth = mean_m_smooth * N, se_M_smooth = se_m_smooth * N)  # count-scale, see header note

b_levels <- sort(unique(grid_w$b))
# Use distinct colours for the six displayed feedback levels.
grid_pal <- c(
  "0"    = gray_off_dark,  # feedback off
  "0.5"  = col_low,        # locked main-text value
  "0.75" = col_on,
  "0.9"  = "#8E6FAE",
  "1"    = "#D4A017",
  "1.3"  = col_high
)
grid_pal <- grid_pal[as.character(b_levels)]  # order/subset to match whatever b's actually ended up in the data
# Use short direct labels and fuller wording in the legend.
grid_labels_plot <- setNames(sprintf("b = %.2f", b_levels), b_levels)
grid_labels_plot[as.character(0)] <- "b = 0"

grid_labels_legend <- setNames(sprintf("b = %.2f", b_levels), b_levels)
grid_labels_legend[as.character(0)] <- "Feedback off (b = 0)"

# Show individual chains for b = 0 and b = 0.50.
# Other feedback levels show the mean and interval.
set.seed(3)
n_sample_chains_headline <- 30L
set.seed(11)
n_sample_chains_grid <- 10L

sim3c_raw_grid <- sim3c_raw |>
  filter(b %in% display_b_grid, time_since_shock >= plot_window_grid[1], time_since_shock <= plot_window_grid[2])

headline_chains <- sim3c_raw_grid |>
  filter(b %in% c(0, 0.5)) |>
  group_by(b) |>
  filter(chain %in% sample(unique(chain), n_sample_chains_headline)) |>
  ungroup()

other_grid_chains <- sim3c_raw_grid |>
  filter(!b %in% c(0, 0.5)) |>
  group_by(b) |>
  filter(chain %in% sample(unique(chain), n_sample_chains_grid)) |>
  ungroup()

grid_chains <- bind_rows(headline_chains, other_grid_chains) |>
  arrange(b, chain, time_since_shock) |>
  group_by(b, chain) |>
  mutate(m = roll_mean(m)) |>
  ungroup() |>
  filter(!is.na(m)) |>
  mutate(M = m * N)  # count-scale version for plotting, see grid_w note above

# grid_w includes b=0 again, so the pre-shock baseline is computed
# directly from those rows -- feedback doesn't change pre-shock dynamics,
# so b=0 is as valid a source for this reference line as any other b.
m_ref_grid <- mean(grid_w$mean_m_smooth[grid_w$b == 0 & grid_w$time_since_shock < 0], na.rm = TRUE)
M_ref_grid <- m_ref_grid * N
x_max_grid <- plot_window_grid[2]

grid_end <- grid_w |> filter(!is.na(mean_M_smooth)) |> group_by(b) |>
  filter(time_since_shock == max(time_since_shock)) |> ungroup() |>
  arrange(mean_M_smooth)
# Spread out end-labels that would otherwise collide (several b's plateau
# close together). label_gap is the original m-scale value (0.045) x N=9,
# now that grid_end is on the count scale.
label_gap <- 0.045 * N
if (nrow(grid_end) > 1) {
  y <- grid_end$mean_M_smooth
  for (i in 2:length(y)) {
    if (y[i] - y[i - 1] < label_gap) y[i] <- y[i - 1] + label_gap
  }
  grid_end$label_y <- y
} else {
  grid_end$label_y <- grid_end$mean_M_smooth
}

# Local tweaks, this panel only: mean lines thinned slightly from the
# shared main_line_width (6 lines sharing the panel reads better a touch
# thinner); end-of-line "b = ..." labels enlarged for legibility.
grid_line_width <- main_line_width * 0.82
grid_label_size <- 4.4

# Include both mean intervals and individual chains in the plotted range.
grid_y_lo <- min(
  grid_w$mean_M_smooth - interval_mult * grid_w$se_M_smooth,
  grid_chains$M,
  na.rm = TRUE
)
grid_y_hi <- max(
  grid_w$mean_M_smooth + interval_mult * grid_w$se_M_smooth,
  grid_chains$M,
  na.rm = TRUE
)
grid_y_pad <- 0.05 * (grid_y_hi - grid_y_lo)
grid_ylim <- c(grid_y_lo - grid_y_pad, grid_y_hi + grid_y_pad)

pGrid <- ggplot(grid_w, aes(x = time_since_shock, y = mean_M_smooth, colour = factor(b), group = factor(b))) +
  geom_hline(yintercept = M_ref_grid, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = grid_chains,
    aes(x = time_since_shock, y = M, group = interaction(chain, b), colour = factor(b)),
    inherit.aes = FALSE, alpha = spaghetti_alpha, linewidth = spaghetti_width
  ) +
  geom_ribbon(
    aes(ymin = mean_M_smooth - interval_mult * se_M_smooth, ymax = mean_M_smooth + interval_mult * se_M_smooth, fill = factor(b)),
    alpha = ribbon_alpha, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = grid_line_width) +
  geom_text(
    data = grid_end, aes(x = x_max_grid + 30, y = label_y, label = grid_labels_plot[as.character(b)], colour = factor(b)),
    hjust = 0, size = grid_label_size, fontface = "plain", show.legend = FALSE
  ) +
  coord_cartesian(xlim = c(plot_window_grid[1], x_max_grid), ylim = grid_ylim, clip = "off") +
  scale_colour_manual(values = grid_pal, labels = grid_labels_legend, name = NULL) +
  scale_fill_manual(values = grid_pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.11))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +  # 6 entries now: wrap to 2 rows so labels don't crowd
  labs(
    title = panel_title("A", "Feedback strength changes whether the system recovers"),
    x = "Steps since perturbation", y = "Mean number of active symptoms"
  ) +
  theme_pub(base_size = base_sz) +
  theme(
    # Right margin bumped up: end-of-line labels are considerably larger
    # now and need headroom to avoid clipping, especially the long
    # "Feedback off (b = 0)" label -- check the rendered PNG for clipping.
    plot.margin = margin(5.5, 88, 5.5, 5.5),
    legend.position = "bottom",
    legend.key.width = unit(0.9, "cm"),
    legend.text = element_text(size = legend_text_sz),
    plot.title = element_text(size = panel_title_sz),
    axis.title = element_text(size = axis_title_sz),
    axis.text = element_text(size = axis_text_sz)
  )

sep <- read_csv("res/revision_2026/supp_history/history_separation_summary.csv", show_col_types = FALSE)

# gray_off_dark in place of col_off, same fix as panel A -- "reconvergent"
# points were rendering in the same washed-out gray that was hard to
# distinguish from the panel background.
regime_pal <- c(convergent = gray_off_dark, `history-dependent` = col_high, `runaway/saturation` = "black")

# Locate the main comparison value b = 0.50 by nearest value, not exact
# equality -- the grid is generated via seq(0, 1.5, by = 0.1), and exact
# `== 0.5` against a seq()-built float risks silently matching nothing due
# to floating-point drift.
locked_b_row <- sep[which.min(abs(sep$b - 0.5)), ]

pB <- ggplot(sep, aes(x = b, y = separation_m)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 0.15, linetype = "dotted", colour = "grey55", linewidth = 0.4) +
  # Connect adjacent grid values while retaining the point-level classification.
  geom_line(colour = "grey45", linewidth = 0.6) +
  geom_pointrange(
    aes(ymin = separation_m - interval_mult * se_separation_m, ymax = separation_m + interval_mult * se_separation_m, colour = regime_flag),
    linewidth = 0.75, size = 0.5
  ) +
  annotate("point", x = locked_b_row$b, y = locked_b_row$separation_m, shape = 21, size = 5.5, colour = col_on, stroke = 1.0) +
  annotate("text", x = locked_b_row$b, y = locked_b_row$separation_m + 0.07, label = "main simulation\nvalue", size = 3.8, colour = col_on, hjust = 0.5, vjust = 0, lineheight = 0.85) +
  # Use clearer display labels without changing the saved regime classification.
  scale_colour_manual(
    values = regime_pal, name = NULL,
    labels = c(
      convergent = "reconvergent",
      `history-dependent` = "initial-state-dependent",
      `runaway/saturation` = "runaway/saturation"
    )
  ) +
  labs(
    title = panel_title("B", "Initial-state dependence emerges at higher feedback"),
    x = "Feedback strength (b)", y = "Difference in late mean symptom activation"
  ) +
  theme_pub(base_size = base_sz) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = legend_text_sz),
    plot.title = element_text(size = panel_title_sz),
    axis.title = element_text(size = axis_title_sz),
    axis.text = element_text(size = axis_text_sz)
  )

# Combine the recovery and initial-state panels.
# Keep separate legends because their colours represent different quantities.
fig3 <- (pGrid / pB) +
  plot_layout(heights = c(1.7, 1))

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave("figs/revision_2026/Figure3_regimes.pdf", fig3, width =12, height = 10)
ggsave("figs/revision_2026/Figure3_regimes.png", fig3, width = 7.8, height = 8.4, dpi = 300)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure3_regimes.pdf (+ .png)\n")
cat("\nThis REPLACES fig3_recovery_feedback.R / FigureS_history_dependence.R\n")
cat("as main-text Figure 3. Those two scripts are kept (not deleted) but are\n")
cat("now superseded for main-text purposes.\n")
