// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "utils.h"
using namespace Rcpp;

void validate_positive_int(int value, const char* name) {
  if (value < 1) {
    Rcpp::stop("%s must be a positive integer", name);
  }
}

void validate_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("%s must contain only finite values", name);
  }
}

void validate_q_vec(const arma::vec& q_vec, int m) {
  if (q_vec.n_elem < static_cast<arma::uword>(m)) {
    Rcpp::stop("q_vec length must be at least m");
  }

  for (int k = 0; k < m; ++k) {
    double q = q_vec(static_cast<arma::uword>(k));
    if (!std::isfinite(q) || q < 1 || q != std::floor(q)) {
      Rcpp::stop("q_vec entries must be positive integers");
    }
  }
}

void validate_constraint_vec(const arma::vec& constraint) {
  if (constraint.n_elem != 3) {
    Rcpp::stop("constraint must have length 3");
  }

  for (arma::uword i = 0; i < constraint.n_elem; ++i) {
    double value = constraint(i);
    if (!std::isfinite(value) || (value != 0.0 && value != 1.0)) {
      Rcpp::stop("constraint entries must be 0/1");
    }
  }
}

void validate_positive_finite_vec(const arma::vec& x,
                                  arma::uword expected_length,
                                  const char* name) {
  if (x.n_elem != expected_length) {
    Rcpp::stop("%s must have length %d", name, static_cast<int>(expected_length));
  }

  for (arma::uword i = 0; i < x.n_elem; ++i) {
    if (!std::isfinite(x(i)) || x(i) <= 0.0) {
      Rcpp::stop("%s entries must be positive finite values", name);
    }
  }
}

double get_positive_finite_slot(Rcpp::S4 obj, const char* slot_name) {
  Rcpp::NumericVector value = obj.slot(slot_name);
  if (value.size() != 1 || !std::isfinite(value[0]) || value[0] <= 0.0) {
    Rcpp::stop("%s must be a positive finite scalar", slot_name);
  }
  return value[0];
}

// [[Rcpp::export]]
arma::mat get_z_matrix_native(const arma::vec& z, int m, int n){

  validate_positive_int(m, "m");
  validate_positive_int(n, "n");
  if (z.n_elem != static_cast<arma::uword>(n)) {
    Rcpp::stop("length of z must equal n");
  }

  arma::mat Zmat = arma::zeros<arma::mat>(m, n);

  for (int j = 0; j < n; j++) {
    double label = z(j);
    if (!std::isfinite(label) || label < 1 || label > m || label != std::floor(label)) {
      Rcpp::stop("cluster labels must be integers in 1:m");
    }
    Zmat(static_cast<arma::uword>(label - 1), j) = 1;
  }
  return(Zmat);
}

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec multivariate_normal_density_native(const arma::mat& x,
                      const arma::rowvec& mean,
                      const arma::mat& sigma,
                      bool logd) {
  int n = x.n_rows;
  int xdim = x.n_cols;

  if (xdim < 1) {
    Rcpp::stop("x must have at least one column");
  }
  if (mean.n_elem != static_cast<arma::uword>(xdim)) {
    Rcpp::stop("mean length must equal the number of columns in x");
  }
  if (sigma.n_rows != sigma.n_cols || sigma.n_rows != static_cast<arma::uword>(xdim)) {
    Rcpp::stop("sigma must be a square matrix with dimension matching x");
  }
  if (!x.is_finite() || !mean.is_finite() || !sigma.is_finite()) {
    Rcpp::stop("x, mean, and sigma must contain only finite values");
  }

  arma::vec out(n);
  arma::mat sigma_chol;
  if (!arma::chol(sigma_chol, sigma)) {
    Rcpp::stop("sigma must be positive definite");
  }

  arma::mat rooti = arma::trans(arma::inv(trimatu(sigma_chol)));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double calculate_ratio_native(double log_denominator, const arma::vec& log_numerator){

  if (!std::isfinite(log_denominator)) {
    Rcpp::stop("log_denominator must be finite");
  }
  if (log_numerator.n_elem < 1) {
    Rcpp::stop("log_numerator must contain at least one value");
  }
  if (!log_numerator.is_finite()) {
    Rcpp::stop("log_numerator must contain only finite values");
  }

  double max_numerator = arma::max(log_numerator);
  double shifted_denominator = log_denominator - max_numerator;

  arma::vec shifted_numerator = log_numerator - max_numerator;
  double ratio = exp(shifted_denominator)/sum(exp(shifted_numerator));

  return(ratio);
}

