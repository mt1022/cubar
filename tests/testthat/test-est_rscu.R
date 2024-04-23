test_that("est_rscu works for internal dataset yeast_cds", {
    rscu <- est_rscu(count_codons(yeast_cds))
    expect_true("data.frame" %in% class(rscu))
    expect_equal(nrow(rscu), 61)
})

test_that('empty codon frequency matrix', {
    expect_equal(nrow(est_rscu(matrix())), 61)
})

test_that('sequence with strong bias', {
    seq <- paste0(rep('AAA', 20), collapse = '')
    seqs <- Biostrings::DNAStringSet(seq)
    rscu <- est_rscu(count_codons(seqs))
    expect_equal(rscu$rscu[rscu$codon == 'AAA'], 1 - 1/22)
})

test_that('sequence with no bias', {
    ct <- get_codon_table()
    ct <- ct[amino_acid != '*']
    seq <- paste0(ct$codon, collapse = '')
    seqs <- Biostrings::DNAStringSet(seq)
    rscu <- est_rscu(count_codons(seqs))
    expect_true(all(rscu$RSCU == 1))
})
