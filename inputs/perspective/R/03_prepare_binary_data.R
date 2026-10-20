# ==============================================================================
# R/03_prepare_binary_data.R
#
# Dichotomize PHQ-9 symptoms for the binary Ising analyses.
#
# Requires:
#   R/00_setup.R
#   R/01_prepare_analysis_data.R
#   R/02_model_helpers.R
# ==============================================================================

threshold <- 1L

analysis_df_bin <- analysis_df

analysis_df_bin[symptoms] <- lapply(
  analysis_df_bin[symptoms],
  to_binary01,
  thr = threshold
)

distinct_counts <- vapply(
  analysis_df_bin[symptoms],
  function(x) length(unique(x[!is.na(x)])),
  integer(1)
)

if (!all(distinct_counts <= 2L)) {
  stop(
    "Binarization failed for: ",
    paste(names(distinct_counts)[distinct_counts > 2L], collapse = ", ")
  )
}
