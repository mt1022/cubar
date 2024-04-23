test_that("works for internal dataset yeast_cds and yeast_half_life", {
  csc <- est_csc(yeast_cds, yeast_half_life)
  expect_equal(sum(csc$csc > 0), 24)
})
