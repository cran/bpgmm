test_that("package citation includes the methodology paper DOI", {
  cites <- suppressWarnings(citation("bpgmm"))
  dois <- vapply(cites, function(x) {
    doi <- x$doi
    if (is.null(doi)) {
      NA_character_
    } else {
      doi
    }
  }, character(1))

  expect_true("10.1007/s00357-021-09391-8" %in% dois)
})
