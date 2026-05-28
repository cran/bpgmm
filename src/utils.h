#ifndef BPGMM_UTILS_H
#define BPGMM_UTILS_H

#include <RcppArmadillo.h>

arma::vec multivariate_normal_density_native(const arma::mat &x,
                      const arma::rowvec &mean,
                      const arma::mat &sigma,
                      bool logd);

arma::mat get_z_matrix_native(const arma::vec &z, int m, int n);

double calculate_ratio_native(double log_denominator, const arma::vec &log_numerator);

void validate_positive_int(int value, const char* name);

void validate_finite_matrix(const arma::mat &x, const char* name);

void validate_q_vec(const arma::vec &q_vec, int m);

void validate_constraint_vec(const arma::vec &constraint);

void validate_positive_finite_vec(const arma::vec &x,
                                  arma::uword expected_length,
                                  const char* name);

double get_positive_finite_slot(Rcpp::S4 obj, const char* slot_name);

Rcpp::IntegerVector update_post_z_native(arma::mat X,
                                 int m,
                                 int n,
                                 Rcpp::S4 theta_y_list);

Rcpp::List update_latent_scores_native(arma::mat X,
                               Rcpp::S4 theta_y_list,
                               arma::vec z,
                               arma::vec clus_ind,
                               arma::vec q_vec);

double evaluate_prior_psi_native(Rcpp::List psy,
                         int p,
                         int m,
                         double delta,
                         double bbeta,
                         arma::vec constraint,
                         arma::vec clus_ind);

double evaluate_prior_lambda_native(int p,
                            int m,
                            double alpha2,
                            arma::vec q_vec,
                            Rcpp::List psy,
                            Rcpp::List lambda,
                            arma::vec constraint,
                            arma::vec clus_ind);

Rcpp::List calculate_cxy_native(int m,
                         int n,
                         Rcpp::S4 hparam,
                         Rcpp::S4 theta_y_list,
                         arma::vec z,
                         arma::vec q_vec,
                         arma::mat X);

Rcpp::List calculate_post_lambda_psi_native(int m,
                                   int p,
                                   Rcpp::S4 hparam,
                                   Rcpp::List cxy_list,
                                   Rcpp::S4 theta_y_list,
                                   arma::vec q_vec,
                                   arma::vec constraint);

Rcpp::S4 update_hyperparameter_native(int m,
                               int p,
                               int q,
                               Rcpp::S4 hparam,
                               Rcpp::S4 theta_y_list,
                               arma::vec d_vec,
                               arma::vec s_vec);

#endif
