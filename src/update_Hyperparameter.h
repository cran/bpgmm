#ifndef __BPGMM_UPDATEHYPERPARAMETER__
#define __BPGMM_UPDATEHYPERPARAMETER__

#include <RcppArmadillo.h>

Rcpp::S4 update_Hyperparameter(
    int m,
    int p,
    int q,
    Rcpp::S4 hparam,
    Rcpp::S4 thetaYList,
    arma::vec dVec,
    arma::vec sVec
  )


#endif
