#ifndef BPGMM_CALCULATE_POST_LAMBDA_PSY_H
#define BPGMM_CALCULATE_POST_LAMBDA_PSY_H

#include <RcppArmadillo.h>
#include "utils.h"

Rcpp::List calculate_post_lambda_psi_native(int m,
                                   int p,
                                   Rcpp::S4 hparam,
                                   Rcpp::List cxy_list,
                                   Rcpp::S4 theta_y_list,
                                   arma::vec q_vec,
                                   arma::vec constraint);


#endif
