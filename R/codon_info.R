#' get codon table by NCBI gene code ID
#'
#' \code{get_codon_table} creates a codon table based on the given id of genetic code in NCBI.
#'
#' @param gcid a string of genetic code id. run `show_codon_tables()` to see available codon tables.
#' @returns a `data.table` with four columns: aa_code, amino_acid, codon, and subfam.
#' @importFrom data.table ':='
#' @export
#'
#' @examples
#' # Standard genetic code
#' get_codon_table()
#'
#' # Vertebrate Mitochondrial genetic code
#' get_codon_table(gcid = '2')
get_codon_table <- function(gcid = '1'){
    amino_acid <- aa_code <- subfam <- codon <- . <- NULL # due to NSE notes in R CMD check
    codon_table <- Biostrings::getGeneticCode(gcid, as.data.frame = TRUE)
    data.table::setDT(codon_table, keep.rownames = 'codon')
    data.table::setnames(codon_table, 'AA', 'aa_code')
    codon_table[, amino_acid := Biostrings::AMINO_ACID_CODE[aa_code]]
    codon_table[aa_code == '*', amino_acid := '*']
    codon_table[, subfam := paste(amino_acid, substr(codon, 1, 2), sep = '_')]
    return(codon_table[, .(aa_code, amino_acid, codon, subfam)])
}


#' create custom codon table from a data frame
#'
#' \code{create_codon_table} creates codon table from data frame of aa to codon mapping.
#'
#' @param aa2codon a data frame with two columns: amino_acid (Ala, Arg, etc.) and codon.
#' @returns a `data.table` with four columns: aa_code, amino_acid, codon, and subfam.
#' @importFrom data.table ':='
#' @export
#'
#' @examples
#' head(aa2codon)
#' create_codon_table(aa2codon = aa2codon)
create_codon_table <- function(aa2codon){
    aa_code <- amino_acid <- subfam <- codon <- . <- NULL # due to NSE notes in R CMD check
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
#'
#' \code{show_codon_tables} print a table of available genetic code from NCBI through
#' `Biostrings::GENETIC_CODE_TABLE`.
#' @returns NULL
#' @export
#'
#' @examples
#' # print available NCBI codon table IDs and descriptions.
#' show_codon_tables()
show_codon_tables <- function(){
    cat(sprintf('%2s: %s',
                Biostrings::GENETIC_CODE_TABLE$id,
                Biostrings::GENETIC_CODE_TABLE$name), sep = '\n')
}
