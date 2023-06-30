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


#' create codon table from data frame of aa to codon mapping
create_codon_table <- function(aa2codon){
    codon_table <- data.table::as.data.table(aa2codon)
    aacode <- c(
        `*` = "*", Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
        Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I", Leu = "L",
        Lys = "K", Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T",
        Trp = "W", Tyr = "Y", Val = "V")
    codon_table[, aa_code := unname(aacode[amino_acid])]
    codon_table[, subfam := paste(amino_acid, substr(codon, 1, 2), sep = '_')]
    return(codon_table[, .(aa_code, amino_acid, codon, subfam)])
}


#' show available codon tables
show_codon_tables <- function(){
    cat(sprintf('%2s: %s',
                Biostrings::GENETIC_CODE_TABLE$id,
                Biostrings::GENETIC_CODE_TABLE$name), sep = '\n')
}
