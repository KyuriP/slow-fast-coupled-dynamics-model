#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Each row of states is one complete binary profile. Inverse-CDF sampling
// is exact conditional on each person's P and avoids burn-in assumptions.
// [[Rcpp::export]]
arma::mat sim5_sample(const arma::mat& states, const arma::vec& base,
                     const arma::vec& field, const arma::vec& context,
                     const arma::vec& uniforms) {
  arma::mat out(context.n_elem, states.n_cols);
  for (arma::uword k = 0; k < context.n_elem; ++k) {
    arma::vec eta = base + context(k) * field;
    arma::vec prob = arma::exp(eta - eta.max());
    double target = uniforms(k) * arma::accu(prob);
    double acc = 0.0;
    arma::uword chosen = prob.n_elem - 1;
    for (arma::uword s = 0; s < prob.n_elem; ++s) {
      acc += prob(s);
      if (acc >= target) { chosen = s; break; }
    }
    out.row(k) = states.row(chosen);
  }
  return out;
}

double sim5_loglik(const arma::mat& X, const arma::vec& y,
                  const arma::vec& b) {
  arma::vec eta = X * b;
  double val = 0;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    double z = eta(i);
    val += y(i) * z - (z > 0 ? z + std::log1p(std::exp(-z))
                              : std::log1p(std::exp(z)));
  }
  return val;
}

// Unpenalized logistic regressions, with line search and explicit failure.
// [[Rcpp::export]]
Rcpp::List sim5_fit_cpp(const arma::mat& Y, const arma::vec& context,
                   bool adjust) {
  const int n = Y.n_rows, p = Y.n_cols, d = p + (adjust ? 1 : 0);
  arma::mat beta(p, p, arma::fill::zeros);
  double worst_score = 0;
  double largest_coefficient = 0;
  bool all_ok = true;
  int max_iterations = 0;
  // A zero cell in any symptom-pair table gives a separating direction
  // in a nodewise regression. A small score alone cannot detect this.
  for (int i = 0; i < p; ++i) for (int j = i + 1; j < p; ++j) {
    int counts[4] = {0,0,0,0};
    for (int k = 0; k < n; ++k) counts[int(2*Y(k,i)+Y(k,j))]++;
    for (int c = 0; c < 4; ++c) if (counts[c] == 0) {
      beta.fill(arma::datum::nan);
      return Rcpp::List::create(Rcpp::Named("J")=beta,
        Rcpp::Named("ok")=false, Rcpp::Named("max_score")=NA_REAL,
        Rcpp::Named("iterations")=0,Rcpp::Named("max_coefficient")=NA_REAL,
        Rcpp::Named("reason")="zero symptom-pair cell");
    }
  }
  for (int i = 0; i < p; ++i) {
    arma::mat X(n, d, arma::fill::ones);
    int col = 1;
    for (int j = 0; j < p; ++j) if (j != i) X.col(col++) = Y.col(j);
    if (adjust) X.col(d - 1) = context;
    arma::vec y = Y.col(i), b(d, arma::fill::zeros);
    const double mean_y = arma::mean(y);
    if (mean_y == 0 || mean_y == 1) { all_ok = false; break; }
    b(0) = std::log(mean_y / (1 - mean_y));
    bool converged = false;
    for (int iter = 0; iter < 100; ++iter) {
      max_iterations = std::max(max_iterations, iter + 1);
      arma::vec eta = X * b;
      arma::vec mu = 1.0 / (1.0 + arma::exp(-eta));
      arma::vec score = X.t() * (y - mu);
      double err = arma::abs(score).max();
      if (err < 1e-7) { converged = true; worst_score = std::max(worst_score, err); break; }
      arma::vec weight = mu % (1 - mu);
      arma::mat H = X.t() * (X.each_col() % weight);
      arma::vec step;
      bool solved = arma::solve(step, H, score,
                               arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
      if (!solved || !step.is_finite()) break;
      double ll = sim5_loglik(X, y, b), scale = 1.0;
      arma::vec next = b + step;
      while (sim5_loglik(X, y, next) < ll - 1e-10 && scale > 1e-8) {
        scale *= 0.5; next = b + scale * step;
      }
      b = next;
      if (!b.is_finite() || arma::abs(b).max() > 50) break;
    }
    if (!converged) { all_ok = false; break; }
    largest_coefficient = std::max(largest_coefficient,arma::abs(b).max());
    col = 1;
    for (int j = 0; j < p; ++j) if (j != i) beta(i,j) = b(col++);
  }
  arma::mat J = 0.5 * (beta + beta.t());
  if (!all_ok) J.fill(arma::datum::nan);
  return Rcpp::List::create(Rcpp::Named("J")=J, Rcpp::Named("ok")=all_ok,
                            Rcpp::Named("max_score")=worst_score,
                            Rcpp::Named("iterations")=max_iterations,
                            Rcpp::Named("max_coefficient")=largest_coefficient,
                            Rcpp::Named("reason")=all_ok ? "success" : "numerical failure");
}
