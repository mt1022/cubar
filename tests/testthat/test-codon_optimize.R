test_that("optimization works", {
    optimal_codons <- data.table::data.table(
        aa_code = rep('K', 2), codon = c('AAA', 'AAG'), coef = c(-0.5, 0.5))
    seq <- codon_optimize('AAGAAA', optimal_codons)
    expect_equal(as.character(seq), 'AAAAAA')
})
