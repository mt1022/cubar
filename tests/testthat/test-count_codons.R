test_that("return is matrix", {
    expect_equal(
        dim(count_codons(yeast_cds[1])),
        c(1L, 64L)
    )
    expect_equal(
        dim(count_codons(yeast_cds)),
        c(6600L, 64L)
    )
})


test_that("codon frequency sum is correct", {
    expect_equal(
        rowSums(count_codons(yeast_cds)),
        Biostrings::width(yeast_cds)/3,
        ignore_attr = TRUE
    )
})
