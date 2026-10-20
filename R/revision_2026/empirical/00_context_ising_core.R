# ============================================================
# R/revision_2026/empirical/00_context_ising_core.R
#
# Exact-likelihood Ising models whose parameters depend on person-level
# covariates. This is the engine for the slow-fast model-fitting analyses
# (scripts 01-04).
#
# Model, {0,1} parameterization, p binary symptoms, q = p(p-1)/2 pairs:
#
#   P(X = x | z_n) proportional to exp{ T(x)' theta_n },
#   T(x)      = (x_1..x_p, x_1x_2, ..., x_{p-1}x_p)        length p + q
#   theta_n   = B %*% z_n,  z_n = (1, z_n1, ..., z_nC)     covariate vector
#
# B is a (p + q) x (1 + C) matrix. Column 1 holds baseline activations h and
# interactions J. Column c + 1 holds the effect of covariate c on each
# parameter. Which cells of B are free, fixed at zero, or tied together is
# set by an integer "structure" matrix S of the same shape:
#
#   S[k, c] == 0      -> B[k, c] fixed at 0
#   S[k, c] == m > 0  -> B[k, c] = par[m]   (repeating m ties cells together)
#
# This single parameterization covers every model in the analysis:
#
#   group split, free h and J  : z = (1, High), High moderates h and J
#   group split, common J      : z = (1, High), High moderates h only
#   continuous context, common : z = (1, SRL), one gamma shared by all h_i
#   continuous context, item   : z = (1, SRL), gamma_i per symptom
#   context-moderated coupling : z = (1, SRL), SRL moderates h and J
#
# so all models are fitted and compared on exactly the same likelihood.
#
# Because p = 9 there are 2^9 = 512 states and every partition function is
# computed exactly. With person-specific theta_n each person has their own
# partition function; computations are chunked over persons.
# ============================================================


# ------------------------------------------------------------------------
# State space and features
# ------------------------------------------------------------------------

ci_pair_index <- function(p) {
  which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
}

ci_features <- function(X, pair_index) {
  # X: n x p binary matrix -> n x (p + q) sufficient-statistic matrix
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  cbind(X, X[, pair_index[, 1], drop = FALSE] * X[, pair_index[, 2], drop = FALSE])
}

ci_setup <- function(symptoms) {
  p <- length(symptoms)
  pair_index <- ci_pair_index(p)
  states <- as.matrix(expand.grid(replicate(p, c(0, 1), simplify = FALSE)))
  colnames(states) <- symptoms
  feat_names <- c(
    paste0("h_", symptoms),
    paste0("J_", symptoms[pair_index[, 1]], "_", symptoms[pair_index[, 2]])
  )
  list(
    symptoms = symptoms,
    p = p,
    q = nrow(pair_index),
    pair_index = pair_index,
    states = states,
    F_states = ci_features(states, pair_index),   # 512 x (p + q)
    feat_names = feat_names,
    is_h = c(rep(TRUE, p), rep(FALSE, nrow(pair_index)))
  )
}


# ------------------------------------------------------------------------
# Structure matrices
# ------------------------------------------------------------------------
#
# ci_structure() builds S from a compact description. `covariates` names the
# columns of Z after the intercept. `effects` is a named list, one entry per
# covariate, each one of:
#
#   "h_item"    : separate effect on every h_i
#   "h_common"  : one effect shared by all h_i
#   "hJ_item"   : separate effect on every h_i and every J_ij
#   "J_item"    : separate effect on every J_ij only
#   "none"      : covariate enters the design but has no effect
#
# Baseline h and J (intercept column) are always free.

ci_structure <- function(setup, covariates = character(0), effects = list()) {
  K <- setup$p + setup$q
  S <- matrix(0L, K, 1 + length(covariates),
              dimnames = list(setup$feat_names, c("(Intercept)", covariates)))
  m <- 0L
  S[, 1] <- seq_len(K)
  m <- K
  for (cv in covariates) {
    type <- effects[[cv]]
    if (is.null(type)) type <- "none"
    if (type == "h_item") {
      idx <- which(setup$is_h)
      S[idx, cv] <- m + seq_along(idx); m <- m + length(idx)
    } else if (type == "h_common") {
      S[setup$is_h, cv] <- m + 1L; m <- m + 1L
    } else if (type == "hJ_item") {
      S[, cv] <- m + seq_len(K); m <- m + K
    } else if (type == "J_item") {
      idx <- which(!setup$is_h)
      S[idx, cv] <- m + seq_along(idx); m <- m + length(idx)
    } else if (type != "none") {
      stop("Unknown effect type: ", type)
    }
  }
  attr(S, "npar") <- m
  S
}

ci_par_to_B <- function(par, S) {
  B <- matrix(0, nrow(S), ncol(S), dimnames = dimnames(S))
  on <- S > 0L
  B[on] <- par[S[on]]
  B
}

ci_B_grad_to_par <- function(GB, S) {
  on <- S > 0L
  g <- numeric(attr(S, "npar"))
  agg <- rowsum(GB[on], S[on])
  g[as.integer(rownames(agg))] <- agg[, 1]
  g
}


# ------------------------------------------------------------------------
# Log-likelihood and gradient (chunked over persons)
# ------------------------------------------------------------------------
#
# T_obs : n x (p + q) observed features (ci_features of the data)
# Z     : n x (1 + C) design, first column all ones
# Returns -logLik and its gradient with respect to par (for optim).

ci_negll_grad <- function(par, S, T_obs, Z, setup, chunk = 4000L,
                          want_grad = TRUE) {
  B <- ci_par_to_B(par, S)
  Fs <- setup$F_states
  n <- nrow(Z)
  ll <- 0
  GB <- matrix(0, nrow(B), ncol(B))
  for (start in seq(1L, n, by = chunk)) {
    idx <- start:min(n, start + chunk - 1L)
    Zc <- Z[idx, , drop = FALSE]
    Theta <- B %*% t(Zc)                      # (p+q) x nc
    eta <- Fs %*% Theta                       # 512 x nc
    mx <- apply(eta, 2L, max)
    E <- exp(sweep(eta, 2L, mx))
    den <- colSums(E)
    logZ <- mx + log(den)
    Tc <- T_obs[idx, , drop = FALSE]
    ll <- ll + sum(rowSums(Tc * t(Theta))) - sum(logZ)
    if (want_grad) {
      P <- sweep(E, 2L, den, "/")             # 512 x nc
      Ex <- t(Fs) %*% P                       # (p+q) x nc
      GB <- GB + (t(Tc) - Ex) %*% Zc
    }
  }
  out <- -ll
  if (want_grad) attr(out, "gradient") <- -ci_B_grad_to_par(GB, S)
  out
}


# ------------------------------------------------------------------------
# Fitting
# ------------------------------------------------------------------------

ci_start <- function(S, T_obs, setup) {
  par <- numeric(attr(S, "npar"))
  prev <- pmin(pmax(colMeans(T_obs[, seq_len(setup$p), drop = FALSE]), .005), .995)
  par[S[seq_len(setup$p), 1]] <- qlogis(prev)
  par
}

# Map from par to vec(B): vec(B) = M %*% par.
ci_par_map <- function(S) {
  on <- which(S > 0L)
  M <- matrix(0, length(S), attr(S, "npar"))
  M[cbind(on, S[on])] <- 1
  M
}

# Exact negative log-likelihood, gradient and Hessian in one pass.
# The Ising log-likelihood is concave in its natural parameters, and theta_n
# is linear in par, so -logLik is convex in par: Newton's method is safe.
#
# Hessian of -logLik with respect to vec(B):
#   sum_n  (z_n z_n') (x) Cov_n[T]
# accumulated per pair of design columns (a, b) as
#   F' diag(sum_n P_n z_na z_nb) F  -  sum_n z_na z_nb E_n E_n'.
ci_negll_grad_hess <- function(par, S, T_obs, Z, setup, chunk = 4000L) {
  B <- ci_par_to_B(par, S)
  Fs <- setup$F_states
  K <- nrow(B); C1 <- ncol(B); n <- nrow(Z)
  ll <- 0
  GB <- matrix(0, K, C1)
  HB <- matrix(0, K * C1, K * C1)
  for (start in seq(1L, n, by = chunk)) {
    idx <- start:min(n, start + chunk - 1L)
    Zc <- Z[idx, , drop = FALSE]
    Theta <- B %*% t(Zc)
    eta <- Fs %*% Theta
    mx <- apply(eta, 2L, max)
    E <- exp(sweep(eta, 2L, mx))
    den <- colSums(E)
    Tc <- T_obs[idx, , drop = FALSE]
    ll <- ll + sum(rowSums(Tc * t(Theta))) - sum(mx + log(den))
    P <- sweep(E, 2L, den, "/")
    Ex <- t(Fs) %*% P                          # K x nc
    GB <- GB + (t(Tc) - Ex) %*% Zc
    for (a in seq_len(C1)) for (b in a:C1) {
      w <- Zc[, a] * Zc[, b]
      blk <- crossprod(Fs, Fs * as.numeric(P %*% w)) -
        (Ex * rep(w, each = K)) %*% t(Ex)
      ra <- (a - 1L) * K + seq_len(K); rb <- (b - 1L) * K + seq_len(K)
      HB[ra, rb] <- HB[ra, rb] + blk
      if (a != b) HB[rb, ra] <- HB[rb, ra] + t(blk)
    }
  }
  M <- ci_par_map(S)
  list(value = -ll,
       gradient = -ci_B_grad_to_par(GB, S),
       hessian = t(M) %*% HB %*% M)
}

ci_fit <- function(S, X, Z, setup, start = NULL, label = NULL,
                   hessian = TRUE, maxit = 100L, tol = 1e-8,
                   T_obs = NULL, verbose = FALSE) {
  if (is.null(T_obs)) T_obs <- ci_features(X, setup$pair_index)
  if (is.null(start)) start <- ci_start(S, T_obs, setup)
  par <- start
  cur <- ci_negll_grad_hess(par, S, T_obs, Z, setup)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    step <- tryCatch(solve(cur$hessian, cur$gradient),
                     error = function(e) qr.solve(cur$hessian + diag(1e-6, length(par)),
                                                  cur$gradient))
    t_step <- 1
    repeat {
      cand <- par - t_step * step
      val <- as.numeric(ci_negll_grad(cand, S, T_obs, Z, setup, want_grad = FALSE))
      if (is.finite(val) && val <= cur$value + 1e-10) break
      t_step <- t_step / 2
      if (t_step < 1e-8) break
    }
    par <- cand
    new <- ci_negll_grad_hess(par, S, T_obs, Z, setup)
    dec <- cur$value - new$value
    cur <- new
    if (verbose) message(sprintf("  iter %d  negLL %.4f  max|g| %.2e",
                                 it, cur$value, max(abs(cur$gradient))))
    if (max(abs(cur$gradient)) < tol * max(1, nrow(Z) / 1000) || abs(dec) < 1e-11) {
      converged <- TRUE; break
    }
  }
  vc <- NULL
  if (hessian) vc <- tryCatch(solve(cur$hessian), error = function(e) NULL)
  n <- nrow(Z)
  k <- attr(S, "npar")
  ll <- -cur$value
  list(
    label = label, S = S, par = par, B = ci_par_to_B(par, S),
    vcov = vc, logLik = ll, npar = k, n = n,
    AIC = 2 * k - 2 * ll, BIC = log(n) * k - 2 * ll,
    convergence = as.integer(!converged), iterations = it,
    max_abs_grad = max(abs(cur$gradient))
  )
}

ci_compare <- function(fits) {
  tab <- do.call(rbind, lapply(fits, function(f) data.frame(
    model = f$label, npar = f$npar, logLik = f$logLik,
    AIC = f$AIC, BIC = f$BIC, convergence = f$convergence,
    max_abs_grad = f$max_abs_grad, stringsAsFactors = FALSE)))
  tab$dAIC <- tab$AIC - min(tab$AIC)
  tab$dBIC <- tab$BIC - min(tab$BIC)
  tab
}

ci_lrt <- function(big, small) {
  stat <- max(0, 2 * (big$logLik - small$logLik))
  df <- big$npar - small$npar
  data.frame(larger = big$label, restricted = small$label,
             LR = stat, df = df,
             p = stats::pchisq(stat, df, lower.tail = FALSE),
             stringsAsFactors = FALSE)
}


# ------------------------------------------------------------------------
# K-fold cross-validated log-likelihood (per person, held-out)
# ------------------------------------------------------------------------

# Held-out log-likelihood for every person. Folds are a random split of the
# rows (same seed -> same folds for every model), so models can be compared
# person by person.
ci_cv_person_loglik <- function(S, X, Z, setup, K = 10L, seed = 1L) {
  set.seed(seed)
  T_all <- ci_features(X, setup$pair_index)
  fold <- sample(rep_len(seq_len(K), nrow(Z)))
  out <- numeric(nrow(Z))
  for (k in seq_len(K)) {
    tr <- fold != k
    f <- ci_fit(S, NULL, Z[tr, , drop = FALSE], setup,
                T_obs = T_all[tr, , drop = FALSE], hessian = FALSE)
    Theta <- f$B %*% t(Z[!tr, , drop = FALSE])
    eta <- setup$F_states %*% Theta
    mx <- apply(eta, 2L, max)
    logZ <- mx + log(colSums(exp(sweep(eta, 2L, mx))))
    out[!tr] <- rowSums(T_all[!tr, , drop = FALSE] * t(Theta)) - logZ
  }
  out
}

ci_cv_loglik <- function(S, X, Z, setup, K = 10L, seed = 1L) {
  mean(ci_cv_person_loglik(S, X, Z, setup, K = K, seed = seed))
}

# Paired comparison of two models' held-out log-likelihoods (per person).
ci_cv_paired <- function(ll_a, ll_b, label_a, label_b) {
  d <- ll_a - ll_b
  se <- stats::sd(d) / sqrt(length(d))
  data.frame(model_a = label_a, model_b = label_b,
             mean_diff_per_person = mean(d), se = se,
             ci_lower = mean(d) - 1.96 * se, ci_upper = mean(d) + 1.96 * se,
             stringsAsFactors = FALSE)
}


# ------------------------------------------------------------------------
# Exact sampling and expected statistics
# ------------------------------------------------------------------------
#
# Each person's state is drawn from their own exact 512-state distribution.
# No MCMC, so there is no burn-in or autocorrelation to worry about.

ci_state_probs <- function(B, Z, setup) {
  eta <- setup$F_states %*% (B %*% t(Z))
  P <- exp(sweep(eta, 2L, apply(eta, 2L, max)))
  sweep(P, 2L, colSums(P), "/")              # 512 x n
}

ci_sample <- function(B, Z, setup, chunk = 4000L) {
  n <- nrow(Z)
  draw <- integer(n)
  for (start in seq(1L, n, by = chunk)) {
    idx <- start:min(n, start + chunk - 1L)
    P <- ci_state_probs(B, Z[idx, , drop = FALSE], setup)
    C <- apply(P, 2L, cumsum)
    u <- stats::runif(length(idx))
    draw[idx] <- colSums(sweep(C, 2L, u, "<")) + 1L
  }
  draw <- pmin(draw, nrow(setup$states))
  X <- setup$states[draw, , drop = FALSE]
  rownames(X) <- NULL
  X
}

# Expected sufficient statistics (sum over persons) under the model.
ci_expected_stats <- function(B, Z, setup, chunk = 4000L) {
  n <- nrow(Z)
  tot <- numeric(ncol(setup$F_states))
  for (start in seq(1L, n, by = chunk)) {
    idx <- start:min(n, start + chunk - 1L)
    P <- ci_state_probs(B, Z[idx, , drop = FALSE], setup)
    tot <- tot + rowSums(t(setup$F_states) %*% P)
  }
  tot
}


# ------------------------------------------------------------------------
# Two-group models fitted from sufficient statistics
# ------------------------------------------------------------------------
#
# The four models of the perspective paper (free/equal J x free/equal h) depend
# on the data only through each group's summed features. Fitting them to
# *expected* features gives the population-level (infinite-replication)
# result of the perspective-paper pipeline under any generating model, including the
# expected likelihood-ratio statistics (non-centrality).

ci_two_group_suffstat_fit <- function(stats_list, n_list, setup) {
  p <- setup$p; q <- setup$q; K <- p + q
  Fs <- setup$F_states
  group_ll <- function(theta, s, n) {
    eta <- as.numeric(Fs %*% theta); mx <- max(eta)
    w <- exp(eta - mx); den <- sum(w)
    list(ll = sum(theta * s) - n * (mx + log(den)),
         grad = s - n * as.numeric(crossprod(Fs, w / den)))
  }
  # parameter maps: list of index vectors into par for (h1,J1,h2,J2)
  maps <- list(
    free_J_free_h   = list(h1 = 1:p, J1 = p + 1:q, h2 = K + 1:p, J2 = K + p + 1:q),
    equal_J_free_h  = list(h1 = 1:p, J1 = 2 * p + 1:q, h2 = p + 1:p, J2 = 2 * p + 1:q),
    free_J_equal_h  = list(h1 = 1:p, J1 = p + 1:q, h2 = 1:p, J2 = p + q + 1:q),
    equal_J_equal_h = list(h1 = 1:p, J1 = p + 1:q, h2 = 1:p, J2 = p + 1:q)
  )
  npar <- c(free_J_free_h = 2 * K, equal_J_free_h = 2 * p + q,
            free_J_equal_h = p + 2 * q, equal_J_equal_h = K)
  n_tot <- sum(unlist(n_list))
  fits <- lapply(names(maps), function(m) {
    mp <- maps[[m]]
    start <- numeric(npar[[m]])
    for (g in 1:2) {
      prev <- pmin(pmax(stats_list[[g]][1:p] / n_list[[g]], .005), .995)
      start[mp[[paste0("h", g)]]] <- qlogis(prev)
    }
    fg <- function(par) {
      tot <- 0; gr <- numeric(length(par))
      for (g in 1:2) {
        th <- c(par[mp[[paste0("h", g)]]], par[mp[[paste0("J", g)]]])
        r <- group_ll(th, stats_list[[g]], n_list[[g]])
        tot <- tot + r$ll
        gi <- c(mp[[paste0("h", g)]], mp[[paste0("J", g)]])
        gr <- gr + tabulate_add(gi, r$grad, length(par))
      }
      list(f = -tot, g = -gr)
    }
    opt <- stats::optim(start, function(x) fg(x)$f, function(x) fg(x)$g,
                        method = "BFGS",
                        control = list(maxit = 10000, reltol = 1e-13))
    pr <- lapply(mp, function(ix) opt$par[ix])
    list(model = m, logLik = -opt$value, npar = npar[[m]],
         AIC = 2 * npar[[m]] + 2 * opt$value,
         BIC = log(n_tot) * npar[[m]] + 2 * opt$value,
         h1 = pr$h1, J1 = pr$J1, h2 = pr$h2, J2 = pr$J2,
         convergence = opt$convergence)
  })
  names(fits) <- names(maps)
  fits
}

tabulate_add <- function(index, values, len) {
  out <- numeric(len)
  agg <- rowsum(values, index)
  out[as.integer(rownames(agg))] <- agg[, 1]
  out
}

# Summary of one two-group analysis, in the perspective paper's terms.
ci_two_group_summary <- function(fits, setup) {
  f <- fits
  tab <- data.frame(
    LR_J = 2 * (f$free_J_free_h$logLik - f$equal_J_free_h$logLik),
    LR_h = 2 * (f$free_J_free_h$logLik - f$free_J_equal_h$logLik),
    bic_best = names(f)[which.min(sapply(f, `[[`, "BIC"))],
    aic_best = names(f)[which.min(sapply(f, `[[`, "AIC"))],
    dBIC_commonJ_vs_free = f$equal_J_free_h$BIC - f$free_J_free_h$BIC,
    dBIC_commonJ_vs_commonh = f$equal_J_free_h$BIC - f$free_J_equal_h$BIC,
    sum_dJ_free = sum(f$free_J_free_h$J2 - f$free_J_free_h$J1),
    stringsAsFactors = FALSE
  )
  dh_free <- f$free_J_free_h$h2 - f$free_J_free_h$h1
  dh_eq <- f$equal_J_free_h$h2 - f$equal_J_free_h$h1
  names(dh_free) <- names(dh_eq) <- setup$symptoms
  tab$n_dh_free_positive <- sum(dh_free > 0)
  tab$n_dh_commonJ_positive <- sum(dh_eq > 0)
  tab$dh_slp_free <- dh_free[["slp"]]
  tab$dh_slp_commonJ <- dh_eq[["slp"]]
  tab
}
