test_that("works for internal data yeast_cds", {
  opres <- est_optimal_codons(yeast_cds)
  expect_equal(sum(opres$coef < 0), 27)
})
