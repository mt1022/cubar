test_that("return correct results", {
    seqs <- Biostrings::DNAStringSet(c(
        'CTCCTA',  # 4-fold degenerate sites
        'ATG',     # 1-fold degenerate sites
        'TATTAC',  # 2-fold degenerate sites
        'CTCCTAATGTATTAC'  # combine
    ))
    cf <- count_codons(seqs)
    expect_equal(get_gc4d(cf), c(0.5, NA, NA, 0.5))
})
