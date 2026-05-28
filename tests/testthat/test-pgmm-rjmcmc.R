test_that("pgmm_rjmcmc default constraint is a numeric vector", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_error(
    pgmm_rjmcmc(x, m_init = 1, m_range = c(1, 2), q_new = 1, burn = 0, niter = 1, verbose = FALSE),
    NA
  )
})

test_that("pgmm_rjmcmc accepts paper model labels for constraints", {
  x <- matrix(rnorm(12), nrow = 2)

  fit <- pgmm_rjmcmc(
    x,
    m_init = 1,
    m_range = c(1, 2),
    q_new = 1,
    burn = 0,
    niter = 1,
    constraint = "UUU",
    verbose = FALSE
  )

  expect_equal(fit$constraint_samples[[1]], model_to_constraint("UUU"))
})

test_that("pgmm_rjmcmc honors zero posterior iterations", {
  set.seed(2026)
  x <- cbind(
    matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
    matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
  )

  fit <- pgmm_rjmcmc(
    x,
    m_init = 2,
    m_range = c(1, 3),
    q_new = 1,
    burn = 0,
    niter = 0,
    verbose = FALSE
  )

  expect_length(fit$allocation_samples, 0)
  expect_length(fit$constraint_samples, 0)
  expect_length(fit$alpha1_samples, 0)
})

test_that("pgmm_rjmcmc can run cluster-number and covariance RJMCMC steps", {
  set.seed(7)
  x <- cbind(
    matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
    matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
  )

  fit <- pgmm_rjmcmc(
    x,
    m_init = 2,
    m_range = c(1, 3),
    q_new = 1,
    burn = 0,
    niter = 1,
    m_step = 1,
    v_step = 1,
    verbose = FALSE
  )

  expect_length(fit$allocation_samples, 1)
  expect_length(fit$active_cluster_samples, 1)
  expect_length(fit$constraint_samples, 1)
})

test_that("pgmm_rjmcmc can run split/combine RJMCMC moves", {
  set.seed(1)
  x <- cbind(
    matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
    matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
  )

  fit <- pgmm_rjmcmc(
    x,
    m_init = 2,
    m_range = c(1, 3),
    q_new = 1,
    burn = 0,
    niter = 1,
    m_step = 1,
    split_combine = 1,
    verbose = FALSE
  )

  expect_length(fit$allocation_samples, 1)
  expect_length(fit$active_cluster_samples, 1)
})

test_that("pgmm_rjmcmc validates user-facing inputs", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_error(pgmm_rjmcmc(x, m_init = 0, m_range = c(1, 2), q_new = 1), "m_init")
  expect_error(pgmm_rjmcmc(x, m_init = 1, m_range = c(2, 1), q_new = 1), "m_range")
  expect_error(pgmm_rjmcmc(x, m_init = 1, m_range = c(1, 2), q_new = 0), "q_new")
  expect_error(pgmm_rjmcmc(x, m_init = 1, m_range = c(1, 2), q_new = 1, burn = -1), "burn")
  expect_error(pgmm_rjmcmc(x, m_init = 1, m_range = c(1, 2), q_new = 1, niter = -1), "niter")
  expect_error(pgmm_rjmcmc(x, m_init = 1, m_range = c(1, 2), q_new = 1, constraint = c(0, 0)), "constraint")
  expect_error(pgmm_rjmcmc(x, m_init = 1, m_range = c(1, 2), q_new = 1, m_step = 2), "m_step")
  expect_error(pgmm_rjmcmc(x, m_init = 1, m_range = c(1, 2), q_new = 1, v_step = NA), "v_step")
  expect_error(pgmm_rjmcmc(x, m_init = 1, m_range = c(1, 2), q_new = 1, split_combine = -1), "split_combine")
  expect_error(pgmm_rjmcmc(as.data.frame(x), m_init = 1, m_range = c(1, 2), q_new = 1), "X")
})
