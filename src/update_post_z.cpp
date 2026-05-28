// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include "utils.h"


// [[Rcpp::export]]
Rcpp::IntegerVector update_post_z_native(arma::mat X,
                                 int m,
                                 int n,
                                 Rcpp::S4 theta_y_list) {

  Rcpp::List lambda = theta_y_list.slot("lambda");
  Rcpp::List M = theta_y_list.slot("M");
  Rcpp::List psy = theta_y_list.slot("psy");
  arma::vec tao = theta_y_list.slot("tao");
  arma::uword p = X.n_rows;

  validate_positive_int(m, "m");
  validate_positive_int(n, "n");
  if (X.n_cols != static_cast<arma::uword>(n)) {
    Rcpp::stop("n must equal the number of columns in X");
  }
  if (p < 1) {
    Rcpp::stop("X must have at least one row");
  }
  validate_finite_matrix(X, "X");
  if (tao.n_elem < static_cast<arma::uword>(m)) {
    Rcpp::stop("tao must contain at least m mixture weights");
  }
  if (lambda.size() < m || M.size() < m || psy.size() < m) {
    Rcpp::stop("theta_y_list parameter lists must contain at least m components");
  }

  std::vector<arma::vec> means(m);
  std::vector<arma::mat> root_inverses(m);
  std::vector<double> log_density_constants(m);
  arma::vec log_tao(m);
  const double log2pi = std::log(2.0 * M_PI);

  for (int k = 0; k < m; ++k) {
    double taok = tao(k);
    if (!std::isfinite(taok) || taok <= 0.0) {
      Rcpp::stop("tao must contain positive finite mixture weights");
    }
    log_tao(k) = std::log(taok);

    arma::vec Mk = M(k);
    arma::mat lambdak = lambda(k);
    arma::mat psyk = psy(k);
    if (Mk.n_elem != p) {
      Rcpp::stop("each M component must have length matching nrow(X)");
    }
    if (lambdak.n_rows != p) {
      Rcpp::stop("each lambda component must have nrow matching nrow(X)");
    }
    if (psyk.n_rows != p || psyk.n_cols != p) {
      Rcpp::stop("each psy component must be a square matrix matching nrow(X)");
    }
    if (!Mk.is_finite()) {
      Rcpp::stop("M must contain only finite values");
    }
    validate_finite_matrix(lambdak, "lambda");
    validate_finite_matrix(psyk, "psy");

    arma::mat var = psyk + lambdak * lambdak.t();
    arma::mat sigma_chol;
    if (!arma::chol(sigma_chol, var)) {
      Rcpp::stop("component covariance matrices must be positive definite");
    }

    means[k] = Mk;
    root_inverses[k] = arma::trans(arma::inv(arma::trimatu(sigma_chol)));
    log_density_constants[k] =
      -(static_cast<double>(p) / 2.0) * log2pi +
      arma::sum(arma::log(root_inverses[k].diag()));
  }

  Rcpp::IntegerVector z_one_dim(n);
  arma::vec log_prob(m);
  arma::vec prob(m);

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < m; ++k) {
      arma::vec z = root_inverses[k] * (X.col(i) - means[k]);
      log_prob(k) = log_tao(k) + log_density_constants[k] - 0.5 * arma::dot(z, z);
    }

    double max_log_prob = arma::max(log_prob);
    prob = arma::exp(log_prob - max_log_prob);
    double prob_sum = arma::sum(prob);
    if (!std::isfinite(prob_sum) || prob_sum <= 0.0) {
      Rcpp::stop("allocation probabilities must be positive finite values");
    }
    prob /= prob_sum;

    double draw = R::runif(0.0, 1.0);
    double cumulative = 0.0;
    z_one_dim(i) = m;
    for (int k = 0; k < m; ++k) {
      cumulative += prob(k);
      if (draw <= cumulative) {
        z_one_dim(i) = k + 1;
        break;
      }
    }
  }

  return z_one_dim;
}
