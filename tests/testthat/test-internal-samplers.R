test_theta <- function(n = 4) {
  new(
    "ThetaYList",
    tao = c(0.55, 0.45),
    psy = list(diag(c(1.2, 1.4)), diag(c(1.1, 1.3))),
    M = list(matrix(c(-0.5, 0.25), nrow = 1), matrix(c(0.75, -0.2), nrow = 1)),
    lambda = list(matrix(c(0.2, -0.1), nrow = 2), matrix(c(-0.15, 0.25), nrow = 2)),
    Y = list(
      matrix(c(-0.4, 0.1, 0.3, -0.2), nrow = 1, ncol = n),
      matrix(c(0.2, -0.3, 0.4, 0.1), nrow = 1, ncol = n)
    )
  )
}

test_hparam <- function() {
  new("Hparam", alpha1 = 2, alpha2 = 3, bbeta = 2, delta = 3, ggamma = 1)
}

test_constraints <- list(
  c(TRUE, TRUE, TRUE),
  c(TRUE, TRUE, FALSE),
  c(TRUE, FALSE, TRUE),
  c(TRUE, FALSE, FALSE),
  c(FALSE, TRUE, TRUE),
  c(FALSE, TRUE, FALSE),
  c(FALSE, FALSE, TRUE),
  c(FALSE, FALSE, FALSE)
)

test_that("prior generators produce valid constrained parameters", {
  p <- 2
  m <- 2
  q_vec <- c(1, 1)

  for (constraint in test_constraints) {
    set.seed(100)
    psy <- bpgmm:::generate_prior_psi(p, m, delta = 3, bbeta = 2, constraint)
    expect_length(psy, m)
    expect_true(all(vapply(psy, function(x) all(dim(x) == c(p, p)) && all(diag(x) > 0), logical(1))))

    lambda <- bpgmm:::generate_prior_lambda(p, m, alpha2 = 3, q_vec, psy, constraint)
    expect_length(lambda, m)
    expect_true(all(vapply(lambda, function(x) all(dim(x) == c(p, 1)), logical(1))))
  }
})

test_that("prior evaluators return finite log densities across constraints", {
  p <- 2
  m <- 2
  theta <- test_theta()
  q_vec <- c(1, 1)
  clus_ind <- c(1, 1)

  for (constraint in test_constraints) {
    expect_true(is.finite(bpgmm:::evaluate_prior_psi(theta@psy, p, m, 3, 2, constraint, clus_ind)))
    expect_true(is.finite(bpgmm:::evaluate_prior_lambda(
      p, m, alpha2 = 3, qVec = q_vec, psy = theta@psy, lambda = theta@lambda,
      constraint = constraint, clusInd = clus_ind
    )))
  }
})

test_that("joint prior uses product of allocation probabilities", {
  p <- 2
  m <- 2
  theta <- test_theta()
  theta@tao <- c(0.8, 0.2)
  hparam <- test_hparam()
  q_vec <- c(1, 1)
  clus_ind <- c(1, 1)
  constraint <- c(FALSE, FALSE, FALSE)

  z_same <- c(1, 1, 1, 1)
  z_mixed <- c(1, 1, 2, 2)

  same_prior <- bpgmm:::evaluate_prior(m, p, c(0, 0), hparam, theta, z_same, q_vec, constraint, clus_ind)
  mixed_prior <- bpgmm:::evaluate_prior(m, p, c(0, 0), hparam, theta, z_mixed, q_vec, constraint, clus_ind)

  expect_equal(
    same_prior - mixed_prior,
    sum(log(theta@tao[z_same])) - sum(log(theta@tao[z_mixed])),
    tolerance = 1e-10
  )
})

test_that("proposal generators and evaluators cover all PGMM covariance constraints", {
  p <- 2
  m <- 2
  n <- 4
  q_vec <- c(1, 1)
  x <- matrix(c(-0.7, 0.1, -0.2, 0.4, 0.8, -0.3, 1.1, 0.2), nrow = p)
  z <- c(1, 1, 2, 2)
  theta <- test_theta(n)
  hparam <- test_hparam()
  cxy <- bpgmm:::calculate_cxy(m, n, hparam, theta, z, q_vec, x)

  for (constraint in test_constraints) {
    set.seed(100)
    post <- bpgmm:::calculate_post_lambda_psi(m, p, hparam, cxy, theta, q_vec, constraint)
    expect_named(post, c("lambda", "psy"))
    expect_length(post$lambda, m)
    expect_length(post$psy, m)
    expect_true(all(vapply(post$lambda, function(x) all(dim(x) == c(p, 1)), logical(1))))
    expect_true(all(vapply(post$psy, function(x) all(dim(x) == c(p, p)) && all(diag(x) > 0), logical(1))))

    set.seed(101)
    lambda <- bpgmm:::calculate_proposal_lambda(hparam, theta, cxy, constraint, m, p, q_vec)
    expect_length(lambda, m)
    expect_true(all(vapply(lambda, function(x) all(dim(x) == c(p, 1)), logical(1))))
    expect_true(is.finite(bpgmm:::evaluate_proposal_lambda(hparam, theta, cxy, constraint, lambda, m, q_vec, p)))

    theta_lambda <- theta
    theta_lambda@lambda <- lambda
    cxy_lambda <- bpgmm:::calculate_cxy(m, n, hparam, theta_lambda, z, q_vec, x)

    set.seed(102)
    psy <- bpgmm:::calculate_proposal_psi(hparam, theta_lambda, cxy_lambda, constraint, m, p, q_vec)
    expect_length(psy, m)
    expect_true(all(vapply(psy, function(x) all(dim(x) == c(p, p)) && all(diag(x) > 0), logical(1))))
    expect_true(is.finite(bpgmm:::evaluate_proposal_psi(hparam, theta_lambda, cxy_lambda, constraint, psy, m, p, q_vec)))
  }
})
