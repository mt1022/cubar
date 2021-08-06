codon_table <- function(gcid = '1'){
    codon_table <- Biostrings::getGeneticCode(gcid, as.data.frame = TRUE)
    data.table::setDT(codon_table, keep.rownames = 'codon')
    data.table::setnames(codon_table, 'AA', 'aa_code')
    codon_table[, amino_acid := Biostrings::AMINO_ACID_CODE[Biostrings::GENETIC_CODE]]
    codon_table[aa_code == '*', amino_acid := '*']
    return(codon_table[, .(aa_code, amino_acid, codon)])
}
