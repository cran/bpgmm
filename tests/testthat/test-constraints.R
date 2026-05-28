test_that("constraint helpers map between paper model names and binary codes", {
  expect_equal(model_to_constraint("CCC"), c(1, 1, 1))
  expect_equal(model_to_constraint("CCU"), c(1, 1, 0))
  expect_equal(model_to_constraint("CUC"), c(1, 0, 1))
  expect_equal(model_to_constraint("CUU"), c(1, 0, 0))
  expect_equal(model_to_constraint("UCC"), c(0, 1, 1))
  expect_equal(model_to_constraint("UCU"), c(0, 1, 0))
  expect_equal(model_to_constraint("UUC"), c(0, 0, 1))
  expect_equal(model_to_constraint("UUU"), c(0, 0, 0))

  expect_equal(constraint_to_model(c(1, 1, 1)), "CCC")
  expect_equal(constraint_to_model(c(0, 0, 0)), "UUU")
  expect_equal(bpgmm:::change_constraint_format("(1, 0, 0)"), "CUU")
})

test_that("constraint helpers reject invalid encodings", {
  expect_error(model_to_constraint("ABC"), "one of")
  expect_error(constraint_to_model(c(1, 0)), "length 3")
  expect_error(constraint_to_model(c(1, 0, 2)), "0/1")
})
