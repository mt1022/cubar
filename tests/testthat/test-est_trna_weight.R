test_that("works for internal yeast_trna_gcn", {
  w <- est_trna_weight(yeast_trna_gcn)
  expect_equal(dim(w), c(61, 8))
})

test_that('works for a vector of only two anticodons', {
  w <- est_trna_weight(c('AAA' = 1, 'GAA' = 1))
  expect_equal(w$W[1:2], c(1.2139, 1.5341))
})
