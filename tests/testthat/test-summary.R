test_that("summarize_pgmm_rjmcmc summarizes posterior allocation and models", {
  pgmm_res <- list(
    allocation_samples = list(c(1, 1, 2), c(1, 2, 2), c(1, 1, 2)),
    constraint_samples = list(c(0, 0, 0), c(1, 0, 0), c(1, 0, 0))
  )

  summary <- summarize_pgmm_rjmcmc(pgmm_res, true_cluster = c(1, 1, 2))

  expect_equal(summary$allocation, c(1, 1, 2))
  expect_equal(as.integer(summary$n_clusters["2"]), 3L)
  expect_equal(as.integer(summary$n_constraints["CUU"]), 2L)
  expect_equal(as.integer(summary$n_constraints["UUU"]), 1L)
  expect_equal(summary$ari, 1)
})

test_that("internal allocation summarizer keeps compatibility alias", {
  z_samples <- list(c(1, 1, 2), c(1, 2, 2), c(1, 1, 2))

  expect_equal(bpgmm:::summarize_allocations(z_samples), c(1, 1, 2))
  expect_equal(bpgmm:::summarize_allocations_legacy(z_samples), bpgmm:::summarize_allocations(z_samples))
})

test_that("summarize_pgmm_rjmcmc validates result structure", {
  expect_error(summarize_pgmm_rjmcmc(list()), "allocation_samples")
  expect_error(
    summarize_pgmm_rjmcmc(list(allocation_samples = list(), constraint_samples = list())),
    "allocation_samples"
  )
  expect_error(
    summarize_pgmm_rjmcmc(list(allocation_samples = list(c(1, 2)), constraint_samples = list())),
    "constraint_samples"
  )
})
