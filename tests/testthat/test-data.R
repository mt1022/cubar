test_that("internal data is loaded correctly", {
    # test human_mt is of class DNAStringSet and has a length of 13
    expect_equal(class(human_mt), "DNAStringSet", ignore_attr = TRUE)
    expect_equal(length(human_mt), 13)

    # test yeast_cds is of class DNAStringSet and has a length of 6600
    expect_equal(class(yeast_cds), "DNAStringSet", ignore_attr = TRUE)
    expect_equal(length(yeast_cds), 6600)

    # test that yeast_half_life is of class data.frame with 3888 rows and 3 columns
    expect_true("data.frame" %in% class(yeast_half_life))
    expect_equal(nrow(yeast_half_life), 3888)

    # test that yeast_exp is of class data.frame with 6685 rows and 3 columns
    # and column names are: gene_id, gene_name, fpkm
    expect_true("data.frame" %in% class(yeast_exp))
    expect_equal(nrow(yeast_exp), 6685)
    expect_equal(colnames(yeast_exp), c("gene_id", "gene_name", "fpkm"))

    # test that yeast_trna_gcn is a named vector with 41 elements and
    # all names are valid codons
    expect_true(class(yeast_trna_gcn) %in% c("numeric", "integer", "table"))
    expect_equal(length(yeast_trna_gcn), 41)
    expect_true(all(names(yeast_trna_gcn) %in% c(
        "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC",
        "ATA", "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG",
        "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC",
        "GCA", "GCG", "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
        "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC",
        "TGA", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG",
        "GGT", "GGC", "GGA", "GGG")))

    # test that aa2codon is a data.frame with 64 rows 2 columns: amino_acid, and codon
    # all codons are valid
    # all amino_acids are valid (amino acids are coded like "Ala", "Arg", etc. and stop
    # codons are coded as "*").
    expect_true("data.frame" %in% class(aa2codon))
    expect_equal(nrow(aa2codon), 64)
    expect_equal(ncol(aa2codon), 2)
    expect_true(all(aa2codon$amino_acid %in% c(
        "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
        "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "*")))

})
