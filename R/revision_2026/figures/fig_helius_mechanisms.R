# ============================================================
# R/revision_2026/figures/fig_helius_mechanisms.R
# ============================================================
# Empirical application, mechanism figure: where the observed HELIUS
# statistics fall among HELIUS-sized samples simulated under each fitted
# mechanism (output of empirical/03_sim_mechanisms.R).
#
# Reads:  res/revision_2026/empirical/helius/step2_simulation/
# Writes: figs/revision_2026/Figure_helius_mechanisms.pdf (+ .png)
# ============================================================

source("R/revision_2026/figures/theme_publication.R")

in_dir <- "res/revision_2026/empirical/helius/step2_simulation"
reps <- read.csv(file.path(in_dir, "replications.csv"))
observed <- read.csv(file.path(in_dir, "observed_summary.csv"))

panels <- c(
  LR_J        = "(A) LR test of equal couplings (df = 36)",
  sum_dJ_free = "(B) Summed coupling difference",
  S_signed    = "(C) Global strength difference",
  M           = "(D) Maximum edge difference",
  LR_h        = "(E) LR test of equal activations (df = 9)",
  dh_slp_free = "(F) Sleep activation difference"
)

long <- do.call(rbind, lapply(names(panels), function(v)
  data.frame(mechanism = mechanism_labels[reps$condition],
             panel = panels[[v]], value = reps[[v]])))
long <- long[is.finite(long$value), ]
long$mechanism <- factor(long$mechanism, levels = names(pal_mechanism))
long$panel <- factor(long$panel, levels = panels)
obs <- data.frame(panel = factor(unname(panels), levels = panels),
                  value = unlist(observed[names(panels)]))

g <- ggplot(long, aes(value, fill = mechanism, colour = mechanism)) +
  geom_density(alpha = 0.30, linewidth = 0.3, n = 64) +
  geom_vline(data = obs, aes(xintercept = value), linetype = 2,
             linewidth = secondary_line_width, colour = axis_colour) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_fill_manual(values = pal_mechanism) +
  scale_colour_manual(values = pal_mechanism) +
  labs(x = NULL, y = "Density", fill = NULL, colour = NULL) +
  theme_pub(base_size = 10)

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/revision_2026/Figure_helius_mechanisms.pdf", g, width = 7, height = 7.2)
ggsave("figs/revision_2026/Figure_helius_mechanisms.png", g, width = 7, height = 7.2, dpi = 300)
cat("  figs/revision_2026/Figure_helius_mechanisms.pdf (+ .png)\n")
