#' Calculate ENC
#'
#' \code{get_enc} computes ENC of each CDS
#'
#' @param seqs CDSs, DNAStringSet or object that can be coerced to DNAStringSet
#' @return vector of ENC values, sequence names are used as vector names
get_enc <- function(cf, codon_table, method = 'X12'){
    cf <- count_codons(seqs)
    if(missing(codon_table)){
        codon_table <- get_codon_table(gcid = '1')
    }

    codon_table <- codon_table[!aa_code == '*']
    codon_list <- split(codon_table$codon, codon_table$subfam)

    f_cf <- sapply(codon_list, function(x){
        mx <- cf[, x, drop = FALSE]
        n <- rowSums(mx)
        # p with pseudo count correction
        p <- (mx + 1) / (n + ncol(mx))
        return(rowSums(p^2))
    })
    n_cf <- sapply(codon_list, function(x){
        mx <- cf[, x, drop = FALSE]
        return(rowSums(mx + 1))
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
get_cai <- function(cf, rscu){
    # exclude single codon sub-family
    rscu <- data.table::as.data.table(rscu)
    rscu[, ss := .N, by = .(subfam)]
    rscu <- rscu[ss > 1]

    # codon frequency per CDS
    cf <- cf[, rscu$codon, drop = FALSE]

    # cai
    cai <- exp(cf %*% matrix(log(rscu$w_cai)) / rowSums(cf))
    return(cai[, 1])
}


#' Calculate TAI
#'
#' Calculate tRNA Adaptation Index (TAI) of each CDS
#'
#' TODO test
#'
#' @param seqs CDS sequences
#' @trna_w tRNA weight for each codon, can be generated with `get_trna_weight`.
#' @return a vector of TAI values
get_tai <- function(cf, trna_w){
    # codon frequency per CDS
    cf <- cf[, trna_w$codon, drop = FALSE]

    # cai
    tai <- exp(cf %*% matrix(log(trna_w$w)) / rowSums(cf))
    return(tai[, 1])
}


#' Calculate GC4d
#'
#' Calculate GC content at synonymous position of codons (using four-fold
#' degenerate sites only)
get_gc4d <- function(cf, codon_table){
    if(missing(codon_table)){
        codon_table <- get_codon_table(gcid =)
    }
    codon_table[, ss := .N, by = .(subfam)]
    codon_table <- codon_table[ss == 4]
    codon_table[, gc3 := substr(codon, 3, 3) %in% c('G', 'C')]

    cf <- cf[, codon_table$codon, drop = FALSE]
    n <- rowSums(cf)
    gc <- rowSums(cf[, codon_table[, codon[gc3 == T]], drop = FALSE])
    return(gc/n)
}


#' Calculate Fop
#'
#' Calculate the fraction of optimal codons (Fop) of each CDS
get_fop <- function(seqs, gcid = '1'){
    cf <- count_codons(seqs)
    optimal_codons <- get_optimal_codons(seqs, gcid = gcid)
    op <- optimal_codons[coef > 0 & pvalue < 0.001, codon]
    rowSums(cf[, op]) / rowSums(cf)
}


#' Calculate CSCg
#'
#' Calculate mean CSC of each CDS
get_cscg <- function(seqs, csc){
    cf <- count_codons(seqs)
    cf <- cf[, csc$codon]
    cp <- cf / rowSums(cf)
    cscg <- cp %*% as.matrix(csc$csc)
    setNames(cscg[, 1], rownames(cscg))
}
