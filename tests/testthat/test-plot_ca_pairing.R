test_that("default setting for standard genetic code works", {
  x <- plot_ca_pairing(plot = FALSE)
  expect_equal(dim(x), c(115, 3))
})
