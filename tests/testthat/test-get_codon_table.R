test_that("get_codon_table works", {
    codon_table <- structure(
        list(
            aa_code = c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "*", "*",
                        "C", "C", "*", "W", "L", "L", "L", "L", "P", "P", "P", "P", "H",
                        "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M", "T", "T",
                        "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V", "V", "V",
                        "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G"),
            amino_acid = c("Phe", "Phe", "Leu", "Leu", "Ser", "Ser", "Ser", "Ser", "Tyr",
                           "Tyr", "*", "*", "Cys", "Cys", "*", "Trp", "Leu", "Leu", "Leu",
                           "Leu", "Pro", "Pro", "Pro", "Pro", "His", "His", "Gln", "Gln",
                           "Arg", "Arg", "Arg", "Arg", "Ile", "Ile", "Ile", "Met", "Thr",
                           "Thr", "Thr", "Thr", "Asn", "Asn", "Lys", "Lys", "Ser", "Ser",
                           "Arg", "Arg", "Val", "Val", "Val", "Val", "Ala", "Ala", "Ala",
                           "Ala", "Asp", "Asp", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly"),
            codon = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT",
                      "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC",
                      "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA",
                      "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG",
                      "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT",
                      "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC",
                      "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA",
                      "GGG"),
            subfam = c("Phe_TT", "Phe_TT", "Leu_TT", "Leu_TT", "Ser_TC", "Ser_TC",
                       "Ser_TC", "Ser_TC", "Tyr_TA", "Tyr_TA", "*_TA", "*_TA", "Cys_TG",
                       "Cys_TG", "*_TG", "Trp_TG", "Leu_CT", "Leu_CT", "Leu_CT", "Leu_CT",
                       "Pro_CC", "Pro_CC", "Pro_CC", "Pro_CC", "His_CA", "His_CA", "Gln_CA",
                       "Gln_CA", "Arg_CG", "Arg_CG", "Arg_CG", "Arg_CG", "Ile_AT", "Ile_AT",
                       "Ile_AT", "Met_AT", "Thr_AC", "Thr_AC", "Thr_AC", "Thr_AC", "Asn_AA",
                       "Asn_AA", "Lys_AA", "Lys_AA", "Ser_AG", "Ser_AG", "Arg_AG", "Arg_AG",
                       "Val_GT", "Val_GT", "Val_GT", "Val_GT", "Ala_GC", "Ala_GC", "Ala_GC",
                       "Ala_GC", "Asp_GA", "Asp_GA", "Glu_GA", "Glu_GA", "Gly_GG", "Gly_GG",
                       "Gly_GG", "Gly_GG")),
        row.names = c(NA,-64L),
        class = c("data.table", "data.frame")
    )
    expect_equal(get_codon_table(), codon_table)
})
