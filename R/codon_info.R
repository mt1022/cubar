#' get codon information
get_codon_table <- function(gcid = '1'){
    codon_table <- Biostrings::getGeneticCode(gcid, as.data.frame = TRUE)
    data.table::setDT(codon_table, keep.rownames = 'codon')
    data.table::setnames(codon_table, 'AA', 'aa_code')
    codon_table[, amino_acid := Biostrings::AMINO_ACID_CODE[aa_code]]
    codon_table[aa_code == '*', amino_acid := '*']
    codon_table[, subfam := paste(amino_acid, substr(codon, 1, 2), sep = '_')]
    return(codon_table[, .(aa_code, amino_acid, codon, subfam)])
}

#' show available codon tables
show_codon_tables <- function(){
    cat(sprintf('%2s: %s', GENETIC_CODE_TABLE$id, GENETIC_CODE_TABLE$name), sep = '\n')
}
