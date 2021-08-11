#' compute ENC
#'
#' \code{get_enc} computes ENC of each CDS
#'
#' @param seqs CDSs, DNAStringSet or object that can be coerced to DNAStringSet
#' @return vector of ENC values, sequence names are used as vector names
get_enc <- function(seqs, method = 'X12'){
    seqs <- DNAStringSet(seqs)
    m <- count_codons(seqs)

    codon_info <- get_codon_table()
    codon_info <- codon_info[!aa_code == '*']
    codon_info[, subfam := paste(amino_acid, substr(codon, 1, 2), sep = '_')]
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

