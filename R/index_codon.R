#' Calculate RSCU
#'
#' \code{get_rscu} returns the RSCU value of codons
#'
#' @param seqs CDS sequences
#' @weight a vector of the same length as \code{seqs} that gives different
#' weights to CDSs when count codons. It could be gene expression levels.
#' @pseudo_cnt pseudo count to avoid dividing by zero. This may occur when
#' only a few sequences are available for RSCU calculation.
#' @gcid ID or name of genetic code. Support for non-standard genetic code will
#' be added in the future.
#' @return a data.table
#' @import data.table
#' @references
get_rscu <- function(seqs, weight = 1, pseudo_cnt = 1, gcid = '1'){
    seqs <- Biostrings::DNAStringSet(seqs)
    codon_freq <- colSums(count_codons(seqs) * weight)

    codon_table <- get_codon_table(gcid)
    codon_table <- codon_table[aa_code != '*']
    codon_table[, cts := codon_freq[codon]]
    codon_table[, `:=`(
        rscu = (cts + pseudo_cnt) / sum(cts + pseudo_cnt),
        w_cai = (cts + pseudo_cnt) / max(cts + pseudo_cnt)), by = .(subfam)]
    return(codon_table[])
}
