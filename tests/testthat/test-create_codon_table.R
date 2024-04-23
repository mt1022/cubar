test_that("create_codon_table works", {
    # an example df with amino acid and codon
    aa2codon <- data.frame(amino_acid = c('Ala', 'Ala', 'Ala'),
                           codon = c('TTT', 'TTC', 'AAA'))
    # expected output
    custom_ct <- structure(list(
        aa_code = c("A", "A", "A"),
        amino_acid = c("Ala", "Ala", "Ala"),
        codon = c("TTT", "TTC", "AAA"), subfam = c("Ala_TT",  "Ala_TT", "Ala_AA")),
        row.names = c(NA, -3L), class = c("data.table", "data.frame"))
    expect_equal(create_codon_table(aa2codon), custom_ct)
})
