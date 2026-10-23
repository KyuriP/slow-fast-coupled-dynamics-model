# Simulation 4 network panels for main Figure 5.
# Show the true network, fixed-context estimate, unadjusted estimate and
# context-adjusted estimate with the same node layout.
# Use the single-run example from sim4_raw.rds, with 10,000 participants per arm.
# The metric panels use a separate set of 30 replicates, also with 10,000 per arm.
# Write Figure4_network_trio.pdf/.png. The filename retains its original numbering.

suppressPackageStartupMessages({
  library(qgraph)
})

source("R/revision_2026/figures/theme_publication.R")   # qgraph_pos_col, qgraph_neg_col
source("R/revision_2026/00_parameters_uncentered01.R")   # symptoms, N

raw <- readRDS("res/revision_2026/sim4/sim4_raw.rds")

W_true <- raw$true_omega
W_baseline <- raw$omega_hat$baseline
W_naive <- raw$omega_hat$naive
W_adjusted <- raw$omega_hat$adjusted

# Hide estimated edges with absolute weight below 0.15 in all three estimates.
# All nonzero generating edges exceed this threshold.
plot_min <- 0.15
zero_small <- function(W) { W[abs(W) < plot_min] <- 0; W }

W_true_plot <- W_true
W_baseline_plot <- zero_small(W_baseline)
W_naive_plot <- zero_small(W_naive)
W_adjusted_plot <- zero_small(W_adjusted)

mats <- list(
  True = W_true_plot, Baseline = W_baseline_plot,
  `Symptom-only` = W_naive_plot, `Context-adjusted` = W_adjusted_plot
)
max_edge <- max(sapply(mats, function(m) max(abs(m))))

# Fixed layout from the true network -- same node positions in every
# panel, so a reader can track a given symptom pair across panels.
L <- qgraph(W_true_plot, layout = "spring", DoNotPlot = TRUE)$layout

# Match symptom abbreviations by name so ordering can change.
phq_abbrev <- c(
  anhedonia     = "ANH",
  depressed     = "DEP",
  sleep         = "SLP",
  energy        = "ENE",
  appetite      = "APP",
  guilt         = "GLT",
  concentration = "CON",
  psychomotor   = "MOT",
  suicidal      = "SUI"
)
short_labels <- unname(phq_abbrev[symptoms])
stopifnot(!anyNA(short_labels))  # catch silently if `symptoms` ever changes

# ------------------------------------------------------------------------
# Highlight estimated edges above the plotting threshold when the
# corresponding generating edge is below it. Use the same rule in all panels.
spurious_col <- col_naive

build_edge_colors <- function(W_est, W_true_ref) {
  cmat <- matrix(qgraph_pos_col, nrow(W_est), ncol(W_est))
  cmat[W_est < 0] <- qgraph_neg_col
  is_spurious <- (abs(W_true_ref) < plot_min) & (abs(W_est) >= plot_min)
  cmat[is_spurious] <- spurious_col
  cmat
}

edge_col_baseline <- build_edge_colors(W_baseline_plot, W_true_plot)
edge_col_naive    <- build_edge_colors(W_naive_plot, W_true_plot)
edge_col_adjusted <- build_edge_colors(W_adjusted_plot, W_true_plot)

# Size labels and nodes for the final two-by-two layout.
plot_one <- function(W, title, edge_color_mat = NULL) {
  # Pass the network matrix through qgraph's input argument.
  args <- list(
    input = W, layout = L, maximum = max_edge, fade = TRUE,
    labels = short_labels,
    # Set panel-title size separately from symptom-label size.
    edge.width = 1.3, vsize = 12, label.cex = 1.3,
    title = title, title.cex = 1.7,
    theme = "classic", DoNotPlot = FALSE,
    # Leave enough margin around peripheral nodes.
    mar = c(2, 2.2, 3.2, 2.2)
  )
  if (is.null(edge_color_mat)) {
    args$posCol <- qgraph_pos_col
    args$negCol <- qgraph_neg_col
  } else {
    args$edge.color <- edge_color_mat
  }
  do.call(qgraph, args)
}

# Include the symptom abbreviation key below the network panels.
abbrev_key <- paste(sprintf("%s = %s", phq_abbrev, names(phq_abbrev)), collapse = ", ")

draw_abbrev_strip <- function() {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  # Reduce the key's text size if it exceeds the available width.
  cex_fit <- 1.0
  while (strwidth(abbrev_key, cex = cex_fit) > 0.96 && cex_fit > 0.4) {
    cex_fit <- cex_fit - 0.02
  }
  text(0.5, 0.5, abbrev_key, cex = cex_fit, col = "grey30")
}

# Place the true and fixed-context networks above the unadjusted and adjusted
# networks, so the two estimates from varying context can be compared side by side.
draw_trio <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(3.1, 3.1, 0.34))
  plot_one(W_true_plot, "(A) True coupling")
  plot_one(W_baseline_plot, "(B) Estimated network, fixed P", edge_col_baseline)
  plot_one(W_naive_plot, "(C) Estimated network, P omitted", edge_col_naive)
  plot_one(W_adjusted_plot, "(D) Estimated network, P included", edge_col_adjusted)
  draw_abbrev_strip()
}

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

# Use a near-square canvas for the network panels and abbreviation key.
pdf("figs/revision_2026/Figure4_network_trio.pdf", width = 7.4, height = 7.5)
draw_trio()
dev.off()

png("figs/revision_2026/Figure4_network_trio.png", width = 7.4, height = 7.5, units = "in", res = 220)
draw_trio()
dev.off()

cat(sprintf("True global strength: %.2f | Baseline (fixed P): %.2f | Symptom-only: %.2f | Context-adjusted: %.2f\n",
            sum(abs(W_true[upper.tri(W_true)])), sum(abs(W_baseline[upper.tri(W_baseline)])),
            sum(abs(W_naive[upper.tri(W_naive)])),
            sum(abs(W_adjusted[upper.tri(W_adjusted)]))))

cat("\nDone. Files:\n")
cat("  figs/revision_2026/Figure4_network_trio.pdf (+ .png)\n")

cat("\nNote: symptom labels are PHQ-9-style abbreviations (", paste(short_labels, collapse=", "),
    ").\nOrange edges = estimated couplings on true-zero pairs (absent in the\n",
    "true network) exceeding the plotting threshold in that panel's estimate.\n", sep = "")

cat("\nCheck the rendered PNG: qgraph's edge.color argument is being passed a\n")
cat("full NxN matrix (build_edge_colors()) rather than a vector -- this is a\n")
cat("documented qgraph input format, but if edges render in the wrong color\n")
cat("or qgraph errors on the matrix form, the fallback is to get the edge\n")
cat("list via qgraph(...)$Edgelist and pass edge.color as a vector aligned\n")
cat("to that edge order instead.\n")
