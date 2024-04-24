test_that("return correct results", {
    seqs <- Biostrings::DNAStringSet(c('ATGATT', 'ATC', 'ATA', 'GCGGCA'))
    cf <- count_codons(seqs)
    expect_equal(get_gc3s(cf), c(0, 1, 0, 0.5))
})
