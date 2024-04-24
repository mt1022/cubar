test_that("return correct results", {
  seqs <- Biostrings::DNAStringSet(c('ATGATT', 'ATC', 'ATA', 'GCG'))
  cf <- count_codons(seqs)
  expect_equal(get_gc(cf), c(1/6, 1/3, 0, 1))
})
