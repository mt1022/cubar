#' Retrieve codon table by NCBI genetic code ID
#'
#' \code{get_codon_table} creates a standardized codon table based on genetic 
#' codes cataloged by NCBI. This function provides the mapping between codons 
#' and amino acids for different organisms and organelles, which is essential 
#' for accurate codon usage analysis.
#'
#' @param gcid A character string specifying the NCBI genetic code ID. Use 
#'   \code{show_codon_tables()} to view all available genetic codes and their 
#'   corresponding IDs. Default is "1" (standard genetic code).
#' @return A data.table with four columns:
#'   \itemize{
#'     \item \code{aa_code}: Single-letter amino acid code
#'     \item \code{amino_acid}: Three-letter amino acid abbreviation  
#'     \item \code{codon}: Three-nucleotide codon sequence
#'     \item \code{subfam}: Codon subfamily identifier (amino_acid_XY format)
#'   }
#' @importFrom data.table ':='
#' @export
#' @examples
#' # Standard genetic code (used by most organisms)
#' standard_code <- get_codon_table()
#' head(standard_code)
#'
#' # Vertebrate mitochondrial genetic code
#' mito_code <- get_codon_table(gcid = '2')
#' head(mito_code)
#'
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


#' Create custom codon table from amino acid-codon mapping
#'
#' \code{create_codon_table} generates a codon table from a user-defined data 
#' frame that maps codons to their corresponding amino acids. This function 
#' enables analysis of non-standard or artificial genetic codes not available 
#' in the NCBI genetic code collection.
#'
#' @param aa2codon A data frame with two required columns:
#'   \itemize{
#'     \item \code{amino_acid}: Three-letter amino acid abbreviations (e.g., "Ala", "Arg")
#'     \item \code{codon}: Corresponding three-nucleotide codon sequences
#'   }
#' @return A data.table with four columns:
#'   \itemize{
#'     \item \code{aa_code}: Single-letter amino acid code
#'     \item \code{amino_acid}: Three-letter amino acid abbreviation
#'     \item \code{codon}: Three-nucleotide codon sequence  
#'     \item \code{subfam}: Codon subfamily identifier (amino_acid_XY format)
#'   }
#' @importFrom data.table ':='
#' @export
#' @examples
#' # View the example amino acid to codon mapping
#' head(aa2codon)
#' 
#' # Create a custom codon table
#' custom_table <- create_codon_table(aa2codon = aa2codon)
#' head(custom_table)
#'
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


#' Display available genetic code tables
#'
#' \code{show_codon_tables} displays a formatted list of all genetic code tables 
#' available from NCBI, showing their ID numbers and descriptive names. This 
#' function helps users identify the appropriate genetic code ID to use with 
#' \code{get_codon_table()}.
#'
#' @return No return value (called for side effects). The function prints a 
#'   formatted table of available genetic codes to the console, with each line 
#'   showing the numeric ID and corresponding organism/organelle description.
#' @export
#' @examples
#' # Display all available NCBI genetic code tables
#' show_codon_tables()
#'
show_codon_tables <- function(){
    cat(sprintf('%2s: %s',
                Biostrings::GENETIC_CODE_TABLE$id,
                Biostrings::GENETIC_CODE_TABLE$name), sep = '\n')
}
