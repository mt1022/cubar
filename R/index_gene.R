#' Calculate ENC
#'
#' \code{get_enc} computes ENC of each CDS
#'
#' @param seqs CDSs, DNAStringSet or object that can be coerced to DNAStringSet
#' @return vector of ENC values, sequence names are used as vector names
get_enc <- function(seqs, method = 'X12'){
    seqs <- Biostrings::DNAStringSet(seqs)
    m <- count_codons(seqs)

    codon_info <- get_codon_table()
    codon_info <- codon_info[!aa_code == '*']
    codon_list <- split(codon_info$codon, codon_info$subfam)

    f_cf <- sapply(codon_list, function(x){
        mx <- m[, x, drop = FALSE]
        n <- rowSums(mx)
        # p with pseudo count correction
        p <- (mx + 1) / (n + ncol(mx))
        return(rowSums(p^2))
    })
    n_cf <- sapply(codon_list, function(x){
        mx <- m[, x, drop = FALSE]
        return(rowSums(mx))
    })
    if(length(seqs) == 1){
        f_cf <- t(f_cf)
        n_cf <- t(n_cf)
    }
    ss <- lengths(codon_list)
    N_single <- sum(ss == 1)
    N_double <- sum(ss == 2) * rowSums(n_cf[, ss == 2, drop = F]) /
        rowSums(n_cf[, ss == 2, drop = F] * f_cf[, ss == 2, drop = F])
    N_triple <- sum(ss == 3) * rowSums(n_cf[, ss == 3, drop = F]) /
        rowSums(n_cf[, ss == 3, drop = F] * f_cf[, ss == 3, drop = F])
    N_quad <- sum(ss == 4) * rowSums(n_cf[, ss == 4, drop = F]) /
        rowSums(n_cf[, ss == 4, drop = F] * f_cf[, ss == 4, drop = F])
    Nc <- N_single + N_double + N_triple + N_quad
    return(Nc)
}

#' Calculate CAI
#'
#' Calculate Codon Adaptation Index (CAI) or each CDS
#'
#' @param seqs CDS sequences
#' @rscu rscu table containing CAI weight for each codon. This table could be
#'   generated with `get_rscu` or you can prepare it manually.
#' @return a vector of CAI values
get_cai <- function(seqs, rscu){
    # exclude single codon sub-family
    rscu <- data.table::as.data.table(rscu)
    rscu[, ss := .N, by = .(subfam)]
    rscu <- rscu[ss > 1]

    # codon frequency per CDS
    seqs <- Biostrings::DNAStringSet(seqs)
    codon_freq <- count_codons(seqs)
    codon_freq <- codon_freq[, rscu$codon, drop = FALSE]

    # cai
    cai <- exp(codon_freq %*% matrix(log(rscu$w_cai)) / rowSums(codon_freq))
    return(cai[, 1])
}
