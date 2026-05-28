// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
S4 update_hyperparameter_native(
    int m,
    int p,
    int q,
    Rcpp::S4 hparam,
    Rcpp::S4 theta_y_list,
    arma::vec d_vec,
    arma::vec s_vec
  ){


  validate_positive_int(m, "m");
  validate_positive_int(p, "p");
  validate_positive_int(q, "q");
  validate_positive_finite_vec(d_vec, 3, "d_vec");
  validate_positive_finite_vec(s_vec, 3, "s_vec");

  double delta   = get_positive_finite_slot(hparam, "delta");
  double ggamma  = get_positive_finite_slot(hparam, "ggamma");
  List M         = theta_y_list.slot("M");
  List lambda    = theta_y_list.slot("lambda");
  List psy       = theta_y_list.slot("psy");

  if (M.size() < m || lambda.size() < m || psy.size() < m) {
    Rcpp::stop("theta_y_list slots must each have length at least m");
  }

  for (int k = 0; k < m; ++k) {
    arma::vec Mk = M(k);
    arma::mat lambdak = lambda(k);
    arma::mat psyk = psy(k);

    if (Mk.n_elem != static_cast<arma::uword>(p)) {
      Rcpp::stop("M vectors must have length p");
    }
    if (lambdak.n_rows != static_cast<arma::uword>(p) ||
        lambdak.n_cols != static_cast<arma::uword>(q)) {
      Rcpp::stop("lambda matrices must have p rows and q columns");
    }
    if (psyk.n_rows != static_cast<arma::uword>(p) ||
        psyk.n_cols != static_cast<arma::uword>(p)) {
      Rcpp::stop("psy matrices must be square with dimension p");
    }
    if (!Mk.is_finite()) {
      Rcpp::stop("M must contain only finite values");
    }
    validate_finite_matrix(lambdak, "lambda");
    validate_finite_matrix(psyk, "psy");
  }

  // update alpha1
  double alpha1Rate = 0;
  arma::vec alpha1RateTemp;
  for(int k = 0; k < m; k++){
    arma::vec Mk = M(k);
    arma::mat psyk = psy(k);
    alpha1RateTemp = 0.5 * trans(Mk) * arma::inv(psyk) * Mk;
    alpha1Rate += alpha1RateTemp(0);
  }
  alpha1Rate += s_vec(0);

  double alpha1Shape = m*p/2.0 + d_vec(0);
  double alpha1Scale = 1.0/alpha1Rate;

  // Rcpp::rgamma(n, shape  , scale), scale = 1/rate;
  arma::vec alpha1Vec =  Rcpp::rgamma(1, alpha1Shape, alpha1Scale);

  double alpha2Rate = 0;
  arma::mat alpha2RateMat = arma::zeros(q,q);
  arma::mat alpha2RateTemp;
  for(int k = 0; k < m; k++){

    arma::mat lambdak = lambda(k);
    arma::mat psyk = psy(k);
    alpha2RateTemp = 0.5 * trans(lambdak) * arma::inv(psyk)  * lambdak;
    alpha2RateMat += alpha2RateTemp;
  }

  alpha2Rate += s_vec(1) + sum(alpha2RateMat.diag());

  double alpha2Shape = q*p/2.0 + d_vec(1);
  double alpha2Scale = 1.0/alpha2Rate;

  //Rcpp::rgamma(n, shape  , scale), scale = 1/rate;
  arma::vec alpha2Vec =  Rcpp::rgamma(1, alpha2Shape, alpha2Scale);


  double bbetaRate = 0;
  arma::mat bbetaRateMat = arma::zeros(p,p);
  arma::mat bbetaRateTemp;
  for(int k = 0; k < m; k++){
    arma::mat psyk = psy(k);
    bbetaRateTemp = arma::inv(psyk);
    bbetaRate  += sum(bbetaRateTemp.diag());
  }
  bbetaRate += d_vec(2);

  double bbetaShape = m*p*delta + d_vec(2);
  double bbetaScale = 1.0/bbetaRate;

  arma::vec bbetaVec =  Rcpp::rgamma(1, bbetaShape, bbetaScale);
  // Creating an object of Hparam class
  S4 newhparam("Hparam");

  // Setting values to the slots
  newhparam.slot("alpha1")  = alpha1Vec(0);
  newhparam.slot("alpha2")  = alpha2Vec(0);
  newhparam.slot("delta")   = delta;
  newhparam.slot("ggamma")  = ggamma;
  newhparam.slot("bbeta")   = bbetaVec(0);

  return(newhparam);
}
