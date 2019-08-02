#ifndef __BPGMM_UPDATEPOSTZ__
#define __BPGMM_UPDATEPOSTZ__

#include <RcppArmadillo.h>

Rcpp::IntegerVector updatePost_Z( arma::mat X, int m, int n, Rcpp::S4 thetaYList);


#endif
