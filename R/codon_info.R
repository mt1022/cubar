library(Biostrings)
library(data.table)

codon_table <- function(id_or_name){
    codon_table <- getGeneticCode(id_or_name, as.data.frame = TRUE)
    setDT(codon_table, keep.rownames = 'codon')
    setnames(codon_table, 'AA', 'aa_code')
    codon_table[, amino_acid := AMINO_ACID_CODE[GENETIC_CODE]]
    codon_table[aa_code == '*', amino_acid := '*']
    return(codon_table[, .(aa_code, amino_acid, codon)])
}
