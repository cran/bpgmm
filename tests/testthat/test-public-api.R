test_that("documented public API remains exported", {
  exports <- getNamespaceExports("bpgmm")

  expect_true(all(c(
    "pgmm_rjmcmc",
    "summarize_pgmm_rjmcmc",
    "constraint_to_model",
    "model_to_constraint"
  ) %in% exports))

  expect_false(any(c(
    "pgmmRJMCMC",
    "summarizePgmmRJMCMC",
    "summerizePgmmRJMCMC"
  ) %in% exports))
})
