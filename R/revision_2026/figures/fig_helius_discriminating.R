# ============================================================
# R/revision_2026/figures/fig_helius_discriminating.R
# ============================================================
# Empirical application, discriminating-predictions figure: observed LR
# statistic for equal couplings vs 95% intervals from each fitted model,
# for other splits of SRL and for the age / ethnicity groupings (output of
# empirical/04_discriminating_tests.R).
#
# Reads:  res/revision_2026/empirical/helius/step3_discriminating_tests/
# Writes: figs/revision_2026/Figure_helius_discriminating.pdf (+ .png)
# ============================================================

source("R/revision_2026/figures/theme_publication.R")

in_dir <- "res/revision_2026/empirical/helius/step3_discriminating_tests"
splits <- read.csv(file.path(in_dir, "a_srl_splits.csv"))
groups <- read.csv(file.path(in_dir, "b_out_of_sample_groupings.csv"))

split_labels <- c(
  "median (bottom 1/2 vs top 1/2)"    = "Halves",
  "terciles (bottom 1/3 vs top 1/3)"  = "Thirds",
  "quartiles (bottom 1/4 vs top 1/4)" = "Quarters",
  "quintiles (bottom 1/5 vs top 1/5)" = "Fifths",
  "Younger vs Older"                  = "Younger vs older",
  "Dutch vs Non-Dutch"                = "Dutch vs non-Dutch"
)
both <- rbind(cbind(panel = "(A) Bottom vs top of Slow Risk Load", splits),
              cbind(panel = "(B) Alternative groupings, predicted by the fitted models", groups))
both$comparison <- factor(split_labels[both$split], levels = rev(split_labels))
both$panel <- factor(both$panel, levels = unique(both$panel))

pred <- both[both$source != "observed", ]
pred$mechanism <- factor(mechanism_labels[pred$source], levels = names(pal_mechanism))
obs <- both[both$source == "observed", ]

dodge <- position_dodge(width = 0.6)
g <- ggplot(pred, aes(y = comparison, colour = mechanism)) +
  geom_errorbarh(aes(xmin = LR_J_q025, xmax = LR_J_q975), height = 0.3,
                 linewidth = main_line_width, position = dodge) +
  geom_point(data = obs, aes(x = LR_J_or_noncentrality, y = comparison),
             inherit.aes = FALSE, shape = 4, size = point_size_main, stroke = 1,
             colour = axis_colour) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = pal_mechanism, drop = TRUE) +
  labs(x = "LR statistic for equal couplings (df = 36)", y = NULL,
       colour = "95% interval under") +
  theme_pub(base_size = 10) +
  theme(panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.25),
        panel.grid.major.y = element_blank())

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/revision_2026/Figure_helius_discriminating.pdf", g, width = 7, height = 4.8)
ggsave("figs/revision_2026/Figure_helius_discriminating.png", g, width = 7, height = 4.8, dpi = 300)
cat("  figs/revision_2026/Figure_helius_discriminating.pdf (+ .png)\n")
