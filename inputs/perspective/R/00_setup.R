# ───────────────────────────────────────────────────────────────────────────────
# 00_setup.R
# Packages, seed, global plotting settings
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(ggrepel)
  library(purrr)
  library(broom)
  library(haven)
  library(skimr)
  library(bootnet)
  library(IsingSampler)
  library(qgraph)
  library(ggpubr)
})

set.seed(1)

group_colors <- c(
  "Low Slow Risk"  = "#1B9E77",
  "High Slow Risk" = "#D95F02",
  "Dutch"          = "#7570B3",
  "Non-Dutch"      = "#E7298A",
  "Younger"        = "#66A61E",
  "Older"          = "#E6AB02"
)

theme_paper <- theme_minimal(base_size = 13, base_family = "Palatino") +
  theme(panel.grid.minor = element_blank())

symptom_vars <- c("anh", "dep", "slp", "ene", "app", "glt", "con", "mot", "sui")
symptoms <- symptom_vars

symptom_order <- c("sui", "slp", "mot", "glt", "ene", "dep", "con", "app", "anh")
symptoms_plot <- if (all(symptom_order %in% symptoms)) symptom_order else symptoms
