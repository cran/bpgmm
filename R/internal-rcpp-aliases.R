#' (internal)
#' @noRd
calculate_cxy <- function(m, n, hparam, theta_y_list, z, q_vec, x) {
  calculate_cxy_native(m, n, hparam, theta_y_list, z, q_vec, x)
}

#' (internal)
#' @noRd
calculate_post_lambda_psi <- function(m, p, hparam, cxy_list, theta_y_list, q_vec, constraint) {
  calculate_post_lambda_psi_native(m, p, hparam, cxy_list, theta_y_list, q_vec, constraint)
}

#' (internal)
#' @noRd
update_post_z_cpp <- function(x, m, n, theta_y_list) {
  update_post_z_native(x, m, n, theta_y_list)
}

#' (internal)
#' @noRd
update_latent_scores_cpp <- function(x, theta_y_list, z, clus_ind, q_vec) {
  update_latent_scores_native(x, theta_y_list, z, clus_ind, q_vec)
}

#' (internal)
#' @noRd
update_hyperparameter <- function(m, p, q, hparam, theta_y_list, d_vec, s_vec) {
  update_hyperparameter_native(m, p, q, hparam, theta_y_list, d_vec, s_vec)
}

#' (internal)
#' @noRd
get_z_mat_cpp <- function(z, m, n) {
  get_z_matrix_native(z, m, n)
}

#' (internal)
#' @noRd
calculate_ratio_cpp <- function(log_denominator, log_numerator) {
  calculate_ratio_native(log_denominator, log_numerator)
}

#' (internal)
#' @noRd
evaluate_prior_psi_cpp <- function(psy, p, m, delta, bbeta, constraint, clus_ind) {
  evaluate_prior_psi_native(psy, p, m, delta, bbeta, constraint, clus_ind)
}

#' (internal)
#' @noRd
evaluate_prior_lambda_cpp <- function(p, m, alpha2, q_vec, psy, lambda, constraint, clus_ind) {
  evaluate_prior_lambda_native(p, m, alpha2, q_vec, psy, lambda, constraint, clus_ind)
}
