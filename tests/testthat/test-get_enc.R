test_that("ENC = 61 for a single sequence with no bias", {
  cf <- matrix(1, nrow = 1, ncol = 64)
  colnames(cf) <- sort(names(Biostrings::GENETIC_CODE))
  expect_equal(get_enc(cf), 61)

  cf <- matrix(0, nrow = 2, ncol = 64)
  colnames(cf) <- sort(names(Biostrings::GENETIC_CODE))
  expect_equal(get_enc(cf), rep(61, 2))
})

test_that("ENC \approx the number of subfamilies for seq with extreme bias", {
    ct <- get_codon_table()
    ct <- ct[ct$aa_code != "*", ]
    ct2 <- ct[!duplicated(ct$subfam), ]
    cf <- matrix(0, nrow = 1, ncol = 64)
    colnames(cf) <- sort(names(Biostrings::GENETIC_CODE))
    cf[, colnames(cf) %in% ct2$codon] <- 100000
    expect_equal(round(get_enc(cf), 2), 23)
})
