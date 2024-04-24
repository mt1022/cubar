test_that("TAI = 1 if there is no bias in tRNA availability", {
  ct <- get_codon_table()
  ct <- ct[ct$aa_code != '*', ]
  ct$w <- rep(1, nrow(ct))
  cf <- matrix(sample.int(64), nrow = 1, ncol = 64)
  colnames(cf) <- sort(names(Biostrings::GENETIC_CODE))
  tai <- get_tai(cf, ct)
  expect_equal(tai, 1)
})
