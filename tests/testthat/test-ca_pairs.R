test_that("default setting for standard genetic code works", {
  x <- ca_pairs()
  expect_equal(dim(x), c(115, 4))
})
