test_that("pgmm_rjmcmc_chains validates parallel arguments", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_error(
    pgmm_rjmcmc_chains(x, m_init = 1, m_range = c(1, 2), q_new = 1, chains = 0),
    "chains"
  )
  expect_error(
    pgmm_rjmcmc_chains(x, m_init = 1, m_range = c(1, 2), q_new = 1, cores = 0),
    "cores"
  )
  expect_error(
    pgmm_rjmcmc_chains(x, m_init = 1, m_range = c(1, 2), q_new = 1, seed = NA),
    "seed"
  )
})

test_that("pgmm_rjmcmc_chains runs deterministic independent chains", {
  set.seed(2026)
  x <- cbind(
    matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
    matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
  )

  fit_a <- pgmm_rjmcmc_chains(
    x,
    m_init = 2,
    m_range = c(1, 3),
    q_new = 1,
    burn = 0,
    niter = 1,
    chains = 2,
    cores = 1,
    seed = 42,
    verbose = FALSE
  )
  fit_b <- pgmm_rjmcmc_chains(
    x,
    m_init = 2,
    m_range = c(1, 3),
    q_new = 1,
    burn = 0,
    niter = 1,
    chains = 2,
    cores = 1,
    seed = 42,
    verbose = FALSE
  )

  expect_s3_class(fit_a, "bpgmm_rjmcmc_chains")
  expect_length(fit_a, 2)
  expect_named(fit_a, c("chain_1", "chain_2"))
  expect_length(attr(fit_a, "chain_seeds"), 2)
  expect_equal(attr(fit_a, "chain_seeds"), attr(fit_b, "chain_seeds"))
  expect_equal(fit_a[[1]]$allocation_samples, fit_b[[1]]$allocation_samples)
  expect_equal(fit_a[[2]]$allocation_samples, fit_b[[2]]$allocation_samples)
})
