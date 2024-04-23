test_that("default settings works correctly", {
  x <- Biostrings::DNAStringSet(c("ATGATGTAA", "ATGATGTAA", "ATGATG"))
  y <- check_cds(x)
  expect_equal(as.character(y), rep('ATG', 2))
})

test_that("works for a single sequence", {
    x <- Biostrings::DNAStringSet(c("ATGATGTAA"))
    y <- check_cds(x)
    expect_equal(as.character(y), 'ATG')
})

test_that("input has more than sequence but only one satisfy all criteria", {
    x <- Biostrings::DNAStringSet(c("ATGATGTAA", "ATGATG"))
    y <- check_cds(x)
    expect_equal(as.character(y), 'ATG')
})
