# ───────────────────────────────────────────────────────────────────────────────
# 02_model_helpers.R
# Small helper functions shared across model/figure scripts
# ───────────────────────────────────────────────────────────────────────────────

# Convert ordinal 0–3 PHQ items to binary 0/1.
# thr = 1 means any endorsement becomes 1.
to_binary01 <- function(v, thr = 1L) {
  if (is.factor(v)) v <- as.character(v)
  v <- suppressWarnings(as.numeric(v))
  ifelse(is.na(v), NA_real_, ifelse(v >= thr, 1, 0))
}

# Fixed circular layout helper, only needed for network-style plots.
circle_layout <- function(p) {
  th <- seq(0, 2 * pi, length.out = p + 1L)[-(p + 1L)]
  cbind(cos(th), sin(th))
}

# Binary endorsement helper for logistic-regression / symptom-association plots.
make_binary <- function(x) {
  as.integer(x > 0)
}