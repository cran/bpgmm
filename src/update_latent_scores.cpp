// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "utils.h"

namespace {

arma::vec as_finite_vector(const Rcpp::List &values,
                           int index,
                           arma::uword expected_length,
                           const char* name) {
  arma::vec out = Rcpp::as<arma::vec>(values[index]);
  if (out.n_elem != expected_length || !out.is_finite()) {
    Rcpp::stop("%s entries must be finite vectors with dimension matching X", name);
  }
  return out;
}

arma::mat as_finite_matrix(const Rcpp::List &values,
                           int index,
                           arma::uword expected_rows,
                           arma::uword expected_cols,
                           const char* name) {
  arma::mat out = Rcpp::as<arma::mat>(values[index]);
  if (out.n_rows != expected_rows || out.n_cols != expected_cols || !out.is_finite()) {
    Rcpp::stop("%s entries have incompatible dimensions or non-finite values", name);
  }
  return out;
}

arma::vec rmvnorm_chol(const arma::vec &mean, const arma::mat &sigma) {
  arma::mat root;
  if (!arma::chol(root, sigma)) {
    Rcpp::stop("latent score covariance must be positive definite");
  }

  arma::vec z(mean.n_elem);
  for (arma::uword i = 0; i < z.n_elem; ++i) {
    z(i) = R::rnorm(0.0, 1.0);
  }

  return mean + root.t() * z;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List update_latent_scores_native(arma::mat X,
                               Rcpp::S4 theta_y_list,
                               arma::vec z,
                               arma::vec clus_ind,
                               arma::vec q_vec) {
  validate_finite_matrix(X, "X");

  const int n = static_cast<int>(X.n_cols);
  const int p = static_cast<int>(X.n_rows);
  const int m = static_cast<int>(q_vec.n_elem);

  validate_positive_int(n, "n");
  validate_positive_int(p, "p");
  validate_positive_int(m, "m");
  validate_q_vec(q_vec, m);

  if (z.n_elem != static_cast<arma::uword>(n)) {
    Rcpp::stop("length of z must equal the number of observations");
  }
  if (clus_ind.n_elem < static_cast<arma::uword>(m)) {
    Rcpp::stop("clus_ind length must be at least the number of q_vec entries");
  }

  Rcpp::List lambda_list = theta_y_list.slot("lambda");
  Rcpp::List psy_list = theta_y_list.slot("psy");
  Rcpp::List mean_list = theta_y_list.slot("M");
  Rcpp::List old_scores = theta_y_list.slot("Y");

  if (lambda_list.size() < m || psy_list.size() < m ||
      mean_list.size() < m || old_scores.size() < m) {
    Rcpp::stop("theta_y_list slots must contain at least length(q_vec) entries");
  }

  Rcpp::List scores = Rcpp::clone(old_scores);

  for (int k = 0; k < m; ++k) {
    double active = clus_ind(k);
    if (!std::isfinite(active) || (active != 0.0 && active != 1.0)) {
      Rcpp::stop("clus_ind entries must be 0/1");
    }
    if (active == 0.0) {
      continue;
    }

    const int q = static_cast<int>(q_vec(k));
    arma::mat lambda = as_finite_matrix(lambda_list, k, p, q, "lambda");
    arma::mat psy = as_finite_matrix(psy_list, k, p, p, "psy");
    arma::vec mean = as_finite_vector(mean_list, k, p, "M");

    arma::mat covariance = psy + lambda * lambda.t();
    arma::mat solved;
    if (!arma::solve(solved, covariance, lambda)) {
      Rcpp::stop("component covariance must be solvable");
    }
    arma::mat D = solved.t();
    arma::mat sigma = arma::eye<arma::mat>(q, q) - D * lambda;
    sigma = 0.5 * (sigma + sigma.t());

    arma::mat component_scores(q, n);
    for (int i = 0; i < n; ++i) {
      double label = z(i);
      if (!std::isfinite(label) || label < 1 || label > m || label != std::floor(label)) {
        Rcpp::stop("cluster labels must be integers in 1:m");
      }

      arma::vec score_mean = arma::zeros<arma::vec>(q);
      if (static_cast<int>(label) == k + 1) {
        score_mean = D * (X.col(i) - mean);
      }
      component_scores.col(i) = rmvnorm_chol(score_mean, sigma);
    }

    scores[k] = component_scores;
  }

  return scores;
}
