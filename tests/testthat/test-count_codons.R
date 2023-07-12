test_that("return is matrix - single", {
  expect_equal(
      dim(count_codons(yeast_cds[1])),
      c(1L, 64L)
  )
})

test_that("return is matrix - multiple", {
    expect_equal(
        dim(count_codons(yeast_cds)),
        c(6600L, 64L)
    )
})
