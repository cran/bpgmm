#ifndef __BPGMM_CALCULATECXY__
#define __BPGMM_CALCULATECXY__

#include <RcppArmadillo.h>

Rcpp::List Calculate_Cxy(int m,
                         int n,
                         Rcpp::S4 hparam,
                         Rcpp::S4 thetaYList,
                         arma::vec ZOneDim,
                         arma::vec qVec,
                         arma::mat X);

#endif
