# Simulation 4 metric panels, combined with the network diagrams in main Figure 5.
# Panel E shows the signed sum of estimated couplings.
# Panel F shows absolute estimated coupling on truly uncoupled symptom pairs.
# Summaries use 30 replicates of 10,000 participants per arm.
# Write Figure4_metric_strip.pdf/.png. The filename retains its original numbering.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("R/revision_2026/figures/theme_publication.R")   # theme_pub(), panel_title(), pal_arm

by_arm <- read.csv("res/revision_2026/sim4/sim4_replicated_by_arm.csv") |>
  mutate(arm = factor(arm, levels = c("baseline", "naive", "adjusted")))

se <- function(x) sd(x) / sqrt(length(x))

# Set text sizes for the strip at the manuscript's full text width.
title_size_strip <- 12.0
axis_title_strip <- 11.2
axis_text_strip  <- 10.2
label_size_strip <- 3.4   # geom_text size is in mm, not pt

# Matches the locked Sim 4 results text (10 edges, weights 0.20-0.45);
# not recomputed here to keep this script's only input the replicated CSV.
true_total_coupling <- 3.05

arm_labels <- c(
  baseline = "Fixed P",
  naive = "P omitted",
  adjusted = "P included"
)

summary_by_arm <- by_arm |>
  group_by(arm) |>
  summarise(
    coupling_mean = mean(total_coupling_est),
    coupling_se   = se(total_coupling_est),

    spurious_mean = mean(spurious_abs_coupling),
    spurious_se   = se(spurious_abs_coupling),

    .groups = "drop"
  )

# Panel E. Compare signed total coupling with the generating value, 3.05.
# Signed sums allow positive and negative errors on true-zero edges to cancel.
# Axis limits are fixed for these results. Check them if new results are plotted.
pE <- ggplot(
  summary_by_arm,
  aes(x = arm, y = coupling_mean, colour = arm)
) +
  geom_hline(
    yintercept = true_total_coupling,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = 0.58,
    y = true_total_coupling,
    label = "Data-generating value",
    size = 2.9,
    colour = "grey45",
    hjust = 0,
    vjust = -0.6
  ) +
  geom_pointrange(
    aes(
      ymin = coupling_mean - coupling_se,
      ymax = coupling_mean + coupling_se
    ),
    linewidth = 0.85,
    size = 0.4,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = sprintf("%.2f", coupling_mean)),
    nudge_x = 0.10,
    hjust = 0,
    size = label_size_strip,
    colour = "grey20"
  ) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_discrete(labels = arm_labels) +
  labs(
    title = "(E) Estimated total coupling",
    x = NULL,
    y = expression(paste("Total coupling (", Sigma * hat(omega)[ij], ")"))
  ) +
  theme_pub(base_size = 9.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_size_strip, margin = margin(b = 8)),
    axis.title = element_text(size = axis_title_strip),
    axis.text = element_text(size = axis_text_strip),
    axis.text.x = element_text(size = axis_text_strip)
  )

# Panel F. Sum absolute estimates over pairs with zero generating coupling.
# Check the fixed axis limits if the input results change.
pF <- ggplot(
  summary_by_arm,
  aes(x = arm, y = spurious_mean, colour = arm)
) +
  geom_pointrange(
    aes(
      ymin = spurious_mean - spurious_se,
      ymax = spurious_mean + spurious_se
    ),
    linewidth = 0.65,
    size = 0.35,
    alpha = 0.85
  ) +
  geom_text(
    aes(label = sprintf("%.2f", spurious_mean)),
    nudge_x = 0.10,
    hjust = 0,
    size = label_size_strip,
    colour = "grey20"
  ) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_discrete(labels = arm_labels) +
  labs(
    title = "(F) Spurious coupling among uncoupled pairs",
    x = NULL,
    y = expression(paste("Spurious coupling (", Sigma * group("|", hat(omega)[ij], "|"), ")"))
  ) +
  theme_pub(base_size = 9.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_size_strip, margin = margin(b = 8)),
    axis.title = element_text(size = axis_title_strip),
    axis.text = element_text(size = axis_text_strip),
    axis.text.x = element_text(size = axis_text_strip)
  )

fig4_metric_strip <- pE | pF

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "figs/revision_2026/Figure4_metric_strip.pdf",
  fig4_metric_strip,
  width = 8.4,
  height = 2.5
)

ggsave(
  "figs/revision_2026/Figure4_metric_strip.png",
  fig4_metric_strip,
  width = 8.4,
  height = 3.3,
  dpi = 300
)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure4_metric_strip.pdf (+ .png)\n")
cat("\nThis is now meant to be placed at \\textwidth in LaTeX, same as the\n")
cat("network trio -- update the \\includegraphics line for this file from\n")
cat("width=0.82\\textwidth to width=\\textwidth (or drop the width arg if the\n")
cat("trio's own includegraphics doesn't specify one explicitly).\n")
cat("\nCombine with fig4_network_trio.R's output (Figure4_network_trio.pdf,\n")
cat("panels A-D) as Figure 4's full panel set -- assembled in LaTeX/Overleaf,\n")
cat("not composited in R (qgraph base-R graphics + ggplot don't combine\n")
cat("directly via patchwork without ggplotify/wrap_elements wiring).\n")
