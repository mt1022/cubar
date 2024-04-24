test_that("CAI = 1 for a reference with no bias", {
    cf <- matrix(1, nrow = 1, ncol = 64)
    colnames(cf) <- sort(names(Biostrings::GENETIC_CODE))
    rscu <- est_rscu(cf)
    cf2 <- cf
    cf2[] <- sample.int(64, 64)
    cai <- get_cai(cf2, rscu)
    expect_equal(cai, 1)
})
