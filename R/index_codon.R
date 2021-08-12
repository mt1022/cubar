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

show_ac_pairing <- fuction(gcid='1'){
    codon_table <- get_codon_table(gcid)
    ca_pairing <- codon_table[]
}

#' Calculate tRNA w
#'
#' \code{get_trna_weight} compute the tRNA weight per codon for TAI calculation.
#' This weight reflects relative tRNA availability for each codon.
#'
#' @param trna_level, named vector of tRNA level, one value for each anticodon.
#'   vector names are anticodons.
#' @return data.table of tRNA expression information
#' @import data.table
get_trna_weight <- function(tlevel, gcid = '1'){
    codon_table <- get_codon_table(gcid)
    codon_table[, anticodon := as.character(Biostrings::reverseComplement(
        Biostrings::DNAStringSet(codon_table$codon)))]
    codon_table <- codon_table[aa_code != '*']


    codon_table[, ac_level := trna_level[anticodon]]
    codon_table[is.na(ac_level), ac_level := 0]

}
