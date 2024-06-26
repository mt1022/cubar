#' Optimize codons
#'
#' \code{codon_optimize} takes a coding sequence (without stop codon) and replace
#' each codon to the corresponding synonymous optimal codon.
#'
#' @param seq DNAString, or an object that can be coerced to a DNAString
#' @param optimal_codons table optimze codons as generated by \code{est_optimal_codons}
#' @param gcid genetic code id, default is '1' (standard genetic code)
#' @return a DNAString of the optimized coding sequence
#' @export
#' @examples
#' optimal_codons <- est_optimal_codons(yeast_cds)
#' seq <- 'ATGCTACGA'
#' codon_optimize(seq, optimal_codons)
codon_optimize <- function(seq, optimal_codons, gcid='1'){
    . <- aa_code <- codon <- coef <- NULL
    if(!inherits(seq, 'DNAString')){
        seq <- Biostrings::DNAString(seq)
    }
    pep <- Biostrings::translate(seq, genetic.code = Biostrings::getGeneticCode(id_or_name2 = gcid))
    aas <- strsplit(as.character(pep), '')[[1]]
    optimal_codons <- data.table::as.data.table(optimal_codons)
    optimal_codons <- optimal_codons[order(aa_code, coef), .(aa_code, codon)]
    optimal_codons <- unique(optimal_codons, by = 'aa_code')
    a2c <- optimal_codons$codon
    names(a2c) <- optimal_codons$aa_code
    codons <- a2c[aas]
    return(Biostrings::DNAString(paste0(codons, collapse = '')))
}
