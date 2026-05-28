test_that("Rcpp get_z_mat_cpp matches the R implementation", {
  labels <- c(1, 2, 1, 3)

  expect_equal(
    bpgmm:::get_z_mat_cpp(labels, m = 3, n = length(labels)),
    bpgmm:::get_z_mat_r(labels, m = 3, n = length(labels))
  )
})

test_that("native Rcpp wrappers use snake_case names", {
  native_names <- ls(getNamespace("bpgmm"), pattern = "^[a-z].*_native$")

  expect_true(all(grepl("^[a-z][a-z0-9_]*$", native_names)))
  expect_false(any(grepl("[A-Z]", native_names)))
  expect_setequal(
    native_names,
    c(
      "calculate_cxy_native",
      "calculate_post_lambda_psi_native",
      "calculate_ratio_native",
      "evaluate_prior_lambda_native",
      "evaluate_prior_psi_native",
      "get_z_matrix_native",
      "multivariate_normal_density_native",
      "update_hyperparameter_native",
      "update_latent_scores_native",
      "update_post_z_native"
    )
  )
})

test_that("Rcpp get_z_mat_cpp rejects labels outside 1:m", {
  expect_error(
    bpgmm:::get_z_mat_cpp(c(1, 0, 2), m = 2, n = 3),
    "cluster labels must be integers in 1:m"
  )
  expect_error(
    bpgmm:::get_z_mat_cpp(c(1, 3, 2), m = 2, n = 3),
    "cluster labels must be integers in 1:m"
  )
})

test_that("Rcpp multivariate_normal_density_native matches mvtnorm::dmvnorm", {
  x <- rbind(c(0, 0), c(1, -1))
  mean <- c(0.25, -0.5)
  sigma <- matrix(c(2, 0.4, 0.4, 1.5), nrow = 2)

  expect_equal(
    as.numeric(bpgmm:::multivariate_normal_density_native(x, mean, sigma, TRUE)),
    as.numeric(mvtnorm::dmvnorm(x, mean = mean, sigma = sigma, log = TRUE)),
    tolerance = 1e-10
  )
})

test_that("Rcpp multivariate_normal_density_native validates dimensions and covariance", {
  x <- matrix(c(1, 2), nrow = 1)

  expect_error(
    bpgmm:::multivariate_normal_density_native(x, mean = 0, sigma = diag(2), logd = TRUE),
    "mean length"
  )
  expect_error(
    bpgmm:::multivariate_normal_density_native(x, mean = c(0, 0), sigma = diag(3), logd = TRUE),
    "sigma"
  )
  expect_error(
    bpgmm:::multivariate_normal_density_native(
      x,
      mean = c(0, 0),
      sigma = matrix(c(1, 2, 2, 1), nrow = 2),
      logd = TRUE
    ),
    "positive definite"
  )
})

test_that("Rcpp calculate_ratio_cpp is stable for large log probabilities", {
  expect_equal(
    bpgmm:::calculate_ratio_cpp(1000, c(1000, 1001, 1002)),
    exp(1000 - 1002) / sum(exp(c(1000, 1001, 1002) - 1002)),
    tolerance = 1e-12
  )
})

test_that("Rcpp calculate_ratio_cpp validates finite log inputs", {
  expect_error(
    bpgmm:::calculate_ratio_cpp(0, numeric()),
    "log_numerator"
  )
  expect_error(
    bpgmm:::calculate_ratio_cpp(Inf, c(0, 1)),
    "finite"
  )
  expect_error(
    bpgmm:::calculate_ratio_cpp(0, c(0, NA)),
    "finite"
  )
})

test_that("Rcpp prior evaluators match closed-form reference values", {
  p <- 2
  m <- 2
  alpha2 <- 3
  delta <- 4
  bbeta <- 2
  q_vec <- c(1, 1)
  clus_ind <- c(1, 1)
  psy <- list(diag(c(1.2, 1.4)), diag(c(1.1, 1.3)))
  lambda <- list(matrix(c(0.2, -0.1), nrow = 2), matrix(c(-0.15, 0.25), nrow = 2))

  psi_expected <- sum(dgamma(1 / diag(psy[[1]]), delta, rate = bbeta, log = TRUE)) +
    sum(dgamma(1 / diag(psy[[2]]), delta, rate = bbeta, log = TRUE))
  lambda_expected <- mvtnorm::dmvnorm(lambda[[1]][, 1], rep(0, p), 1 / alpha2 * psy[[1]], log = TRUE) +
    mvtnorm::dmvnorm(lambda[[2]][, 1], rep(0, p), 1 / alpha2 * psy[[2]], log = TRUE)

  expect_equal(
    bpgmm:::evaluate_prior_psi_cpp(psy, p, m, delta, bbeta, c(0, 0, 0), clus_ind),
    psi_expected,
    tolerance = 1e-10
  )
  expect_equal(
    bpgmm:::evaluate_prior_lambda_cpp(p, m, alpha2, q_vec, psy, lambda, c(0, 0, 0), clus_ind),
    lambda_expected,
    tolerance = 1e-10
  )
})

test_that("Rcpp update_post_z_cpp uses log mixture weights", {
  n <- 2000
  theta <- new(
    "ThetaYList",
    tao = c(0.99, 0.01),
    psy = list(matrix(1), matrix(1)),
    M = list(0, 0),
    lambda = list(matrix(0), matrix(0)),
    Y = list(matrix(0, nrow = 1, ncol = n), matrix(0, nrow = 1, ncol = n))
  )

  set.seed(22)
  labels <- bpgmm:::update_post_z_cpp(matrix(0, nrow = 1, ncol = n), m = 2, n = n, theta)

  expect_gt(mean(labels == 1), 0.95)
})

test_that("Rcpp update_post_z_cpp validates dimensions and mixture weights", {
  theta <- new(
    "ThetaYList",
    tao = c(1, 0),
    psy = list(matrix(1), matrix(1)),
    M = list(0, 0),
    lambda = list(matrix(0), matrix(0)),
    Y = list(matrix(0, nrow = 1, ncol = 1), matrix(0, nrow = 1, ncol = 1))
  )

  expect_error(
    bpgmm:::update_post_z_cpp(matrix(0, nrow = 1, ncol = 1), m = 2, n = 1, theta),
    "tao"
  )
  expect_error(
    bpgmm:::update_post_z_cpp(matrix(0, nrow = 1, ncol = 2), m = 2, n = 1, theta),
    "n"
  )
})

test_that("Rcpp update_latent_scores_cpp samples latent score matrices", {
  x <- matrix(c(-1, 0, 1, 2, 3, 4), nrow = 2)
  theta <- new(
    "ThetaYList",
    tao = c(0.5, 0.5),
    psy = list(diag(2), diag(c(1.2, 1.4))),
    M = list(c(0, 0), c(1, 1)),
    lambda = list(matrix(c(0.2, -0.1), nrow = 2), matrix(c(0.1, 0.3), nrow = 2)),
    Y = list(matrix(0, nrow = 1, ncol = 3), matrix(0, nrow = 1, ncol = 3))
  )

  set.seed(2026)
  scores <- bpgmm:::update_latent_scores_cpp(
    x = x,
    theta_y_list = theta,
    z = c(1, 2, 1),
    clus_ind = c(1, 1),
    q_vec = c(1, 1)
  )

  expect_type(scores, "list")
  expect_length(scores, 2)
  expect_true(all(vapply(scores, function(y) identical(dim(y), c(1L, 3L)), logical(1))))
  expect_true(all(is.finite(unlist(scores))))
})

test_that("Rcpp update_latent_scores_cpp validates native inputs", {
  theta <- new(
    "ThetaYList",
    tao = c(1),
    psy = list(diag(2)),
    M = list(c(0, 0)),
    lambda = list(matrix(c(0.1, 0.2), nrow = 2)),
    Y = list(matrix(0, nrow = 1, ncol = 2))
  )
  x <- matrix(c(1, 2, 3, 4), nrow = 2)

  expect_error(
    bpgmm:::update_latent_scores_cpp(x, theta, z = c(1, 2), clus_ind = c(1), q_vec = c(1)),
    "cluster labels"
  )
  expect_error(
    bpgmm:::update_latent_scores_cpp(x, theta, z = c(1, 1), clus_ind = c(1), q_vec = c(0)),
    "q_vec"
  )
})

test_that("Rcpp calculate_cxy validates native inputs", {
  theta <- new(
    "ThetaYList",
    tao = c(1),
    psy = list(diag(2)),
    M = list(c(0, 0)),
    lambda = list(matrix(c(0.1, 0.2), nrow = 2)),
    Y = list(matrix(c(0.3, 0.4), nrow = 1, ncol = 2))
  )
  hparam <- new("Hparam", alpha1 = 2, alpha2 = 3, bbeta = 2, delta = 3, ggamma = 1)
  x <- matrix(c(1, 2, 3, 4), nrow = 2)

  expect_error(
    bpgmm:::calculate_cxy(1, 2, hparam, theta, c(1, 1), c(0), x),
    "q_vec entries must be positive"
  )
  x_bad <- x
  x_bad[1, 1] <- Inf
  expect_error(
    bpgmm:::calculate_cxy(1, 2, hparam, theta, c(1, 1), c(1), x_bad),
    "X must contain only finite values"
  )
})

test_that("Rcpp posterior lambda/psi update validates constraints", {
  theta <- new(
    "ThetaYList",
    tao = c(1),
    psy = list(diag(2)),
    M = list(c(0, 0)),
    lambda = list(matrix(c(0.1, 0.2), nrow = 2)),
    Y = list(matrix(c(0.3, 0.4), nrow = 1, ncol = 2))
  )
  hparam <- new("Hparam", alpha1 = 2, alpha2 = 3, bbeta = 2, delta = 3, ggamma = 1)
  x <- matrix(c(1, 2, 3, 4), nrow = 2)
  cxy <- bpgmm:::calculate_cxy(1, 2, hparam, theta, c(1, 1), c(1), x)

  expect_error(
    bpgmm:::calculate_post_lambda_psi(1, 2, hparam, cxy, theta, c(1), c(1, 0)),
    "constraint must have length 3"
  )
  expect_error(
    bpgmm:::calculate_post_lambda_psi(1, 2, hparam, cxy, theta, c(1), c(1, 0, 2)),
    "constraint entries must be 0/1"
  )
})

test_that("Rcpp hyperparameter update validates native hyperparameters", {
  theta <- new(
    "ThetaYList",
    tao = c(1),
    psy = list(diag(2)),
    M = list(c(0, 0)),
    lambda = list(matrix(c(0.1, 0.2), nrow = 2)),
    Y = list(matrix(c(0.3, 0.4), nrow = 1, ncol = 2))
  )
  hparam <- new("Hparam", alpha1 = 2, alpha2 = 3, bbeta = 2, delta = 3, ggamma = 1)

  expect_error(
    bpgmm:::update_hyperparameter(1, 2, 1, hparam, theta, c(1, 1), c(1, 1, 1)),
    "d_vec must have length 3"
  )
  expect_error(
    bpgmm:::update_hyperparameter(1, 2, 1, hparam, theta, c(1, 1, 1), c(1, NA, 1)),
    "s_vec entries must be positive finite values"
  )
})
