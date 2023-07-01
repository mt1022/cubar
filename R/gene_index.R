#' Calculate ENC
#'
#' \code{get_enc} computes ENC of each CDS
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param codon_table codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @return vector of ENC values, sequence names are used as vector names
get_enc <- function(cf, codon_table = get_codon_table()){
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

    if(nrow(cf) == 1){
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
#' \code{get_cai} calculates Codon Adaptation Index (CAI) of each input CDS
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param rscu rscu table containing CAI weight for each codon. This table could be
#'   generated with `est_rscu` or prepared manually.
#' @returns a named vector of CAI values
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
#' \code{get_tai} calculates tRNA Adaptation Index (TAI) of each CDS
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param trna_w tRNA weight for each codon, can be generated with `est_trna_weight()`.
#' @returns a named vector of TAI values
get_tai <- function(cf, trna_w){
    # codon frequency per CDS
    cf <- cf[, trna_w$codon, drop = FALSE]

    # cai
    tai <- exp(cf %*% matrix(log(trna_w$w)) / rowSums(cf))
    return(tai[, 1])
}


#' GC contents at 4-fold degenerate sites
#'
#' Calculate GC content at synonymous position of codons (using four-fold degenerate sites only).
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @returns a named vector of GC4d values.
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


#' Fraction of optimal codons (Fop)
#'
#' \code{get_fop} calculates the fraction of optimal codons (Fop) of each CDS.
#'
#' @param seqs CDS sequences of all protein-coding genes. One for each gene.
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @returns a named vector of fop values.
get_fop <- function(seqs, codon_table = get_codon_table()){
    cf <- count_codons(seqs)
    optimal_codons <- get_optimal_codons(seqs, ctab = codon_table)
    op <- optimal_codons[coef > 0 & pvalue < 0.001, codon]
    rowSums(cf[, op]) / rowSums(cf)
}


#' Mean Codon Stabilization Coefficients
#'
#' \code{get_cscg} calculates Mean Codon Stabilization Coefficients of each CDS.
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param csc table of Codon Stabilization Coefficients as calculated by `est_csc()`.
#' @returns a named vector of cscg values.
get_cscg <- function(cf, csc){
    cf <- cf[, csc$codon]
    cp <- cf / rowSums(cf)
    cscg <- cp %*% as.matrix(csc$csc)
    setNames(cscg[, 1], rownames(cscg))
}
