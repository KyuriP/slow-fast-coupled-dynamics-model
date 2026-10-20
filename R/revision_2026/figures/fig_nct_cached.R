library(ggplot2)
analysis_tag <- "weighted"
summary_tbl <- read.csv("res/NCT_analysis/nct_mc_summary_weighted.csv")
dir.create("figs", showWarnings = FALSE)
base_plot <- function(data, y_var, title, subtitle, y_lab) {
  ggplot(
    data,
    aes(
      x      = factor(n_per_grp),
      y      = .data[[y_var]],
      colour = factor(threshold_gap),
      group  = factor(threshold_gap)
    )
  ) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ topology) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Sample size per group",
      y = y_lab,
      colour = "Threshold gap"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# 8. Plots
# -----------------------------------------------------------------------------
plot_structure <- base_plot(
  summary_tbl,
  "reject_nwinv",
  sprintf("NCT structure invariance rejection rate (%s)", analysis_tag),
  "Groups differ only in threshold. The dashed line marks 5%.",
  "Rejection rate"
)

plot_strength <- base_plot(
  summary_tbl,
  "reject_glstr",
  sprintf("NCT global strength rejection rate (%s)", analysis_tag),
  "Groups differ only in threshold. The dashed line marks 5%.",
  "Rejection rate"
)

plot_edge_fp <- base_plot(
  summary_tbl,
  "fp_edge_05",
  sprintf("Edge-level false-positive rate (%s)", analysis_tag),
  "Proportion of edges flagged at p < .05. The dashed line marks 5%.",
  "False-positive rate"
)

ggsave(
  sprintf("figs/nct_structure_%s.pdf", analysis_tag),
  plot_structure,
  width = 8.5,
  height = 4.8
)

ggsave(
  sprintf("figs/nct_strength_%s.pdf", analysis_tag),
  plot_strength,
  width = 8.5,
  height = 4.8
)

ggsave(
  sprintf("figs/nct_edge_fp_%s.pdf", analysis_tag),
  plot_edge_fp,
  width = 8.5,
  height = 4.8
)

# -----------------------------------------------------------------------------
