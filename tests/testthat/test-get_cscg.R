test_that("CSCg = 0 if all CSC are 0s", {
    ct <- get_codon_table()
    ct <- ct[ct$aa_code != "*", ]
    ct$csc <- rep(0, nrow(ct))
    cf <- matrix(sample.int(64), nrow = 1, ncol = 64)
    colnames(cf) <- sort(names(Biostrings::GENETIC_CODE))
    # works for a single sequence
    expect_equal(get_cscg(cf, ct), 0)
    # works for multiple sequences
    expect_equal(get_cscg(cf[c(1, 1), ], ct), c(0, 0))
})
