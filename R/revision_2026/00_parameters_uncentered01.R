# Parameters for the nine-symptom, uncentered 0/1 model.
# These are illustrative values, not estimates from participant data.
# The same parameters are used across the revised simulations.

symptoms <- c(
  "anhedonia",
  "depressed",
  "sleep",
  "energy",
  "appetite",
  "guilt",
  "concentration",
  "psychomotor",
  "suicidal"
)

N <- length(symptoms)

tau <- c(
  anhedonia     = -2.40,
  depressed     = -2.30,
  sleep         = -2.00,
  energy        = -1.90,
  appetite      = -2.50,
  guilt         = -2.80,
  concentration = -2.30,
  psychomotor   = -3.00,
  suicidal      = -4.00
) + 1.3  # Fixed shift selected during calibration to avoid floor activation.

gamma <- c(
  anhedonia     = 0.90,
  depressed     = 1.00,
  sleep         = 0.70,
  energy        = 0.80,
  appetite      = 0.70,
  guilt         = 1.00,
  concentration = 0.80,
  psychomotor   = 0.70,
  suicidal      = 0.60
)

omega <- matrix(0, nrow = N, ncol = N)
rownames(omega) <- symptoms
colnames(omega) <- symptoms

add_edge <- function(mat, a, b, w) {
  mat[a, b] <- w
  mat[b, a] <- w
  mat
}

omega <- add_edge(omega, "anhedonia", "depressed", 0.45)
omega <- add_edge(omega, "depressed", "guilt", 0.35)
omega <- add_edge(omega, "depressed", "suicidal", 0.30)
omega <- add_edge(omega, "sleep", "energy", 0.40)
omega <- add_edge(omega, "energy", "concentration", 0.35)
omega <- add_edge(omega, "concentration", "psychomotor", 0.25)
omega <- add_edge(omega, "appetite", "energy", 0.25)
omega <- add_edge(omega, "guilt", "suicidal", 0.25)
omega <- add_edge(omega, "anhedonia", "energy", 0.25)
omega <- add_edge(omega, "anhedonia", "concentration", 0.20)

params_uncentered01 <- list(
  symptoms = symptoms,
  N = N,
  tau = tau,
  gamma = gamma,
  omega = omega
)

dir.create("res/revision_2026", recursive = TRUE, showWarnings = FALSE)
saveRDS(params_uncentered01, "res/revision_2026/params_uncentered01.rds")
