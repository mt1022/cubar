#' Calculate ENC
#'
#' \code{get_enc} computes ENC of each CDS
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param codon_table codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @param level "subfam" (default) or "amino_acid". For which level to determine ENC.
#' @return vector of ENC values, sequence names are used as vector names
#' @export
#' @references
#' - Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87:23-29.
#' - Sun X, Yang Q, Xia X. 2013. An improved implementation of effective number of codons (NC).
#'   Mol Biol Evol 30:191-196.
#' @examples
#' # estimate ENC of yeast genes
#' cf_all <- count_codons(yeast_cds)
#' enc <- get_enc(cf_all)
#' head(enc)
#' hist(enc)
get_enc <- function(cf, codon_table = get_codon_table(), level = 'subfam'){
    aa_code <- NULL # due to NSE notes in R CMD check
    if(!level %in% c('amino_acid', 'subfam')){
      stop('Possible values for `level` are "amino_acid" and "subfam"')
    }
    codon_table <- data.table::as.data.table(codon_table)
    codon_table <- codon_table[!aa_code == '*']
    codon_list <- split(codon_table$codon, codon_table[[level]])

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

    Nc <- N_single <- sum(ss == 1)
    if(sum(ss == 2) > 0){
        N_double <- sum(ss == 2) * rowSums(n_cf[, ss == 2, drop = F]) /
            rowSums(n_cf[, ss == 2, drop = F] * f_cf[, ss == 2, drop = F])
        Nc <- Nc + N_double
    }
    if(sum(ss == 3) > 0){
        N_triple <- sum(ss == 3) * rowSums(n_cf[, ss == 3, drop = F]) /
            rowSums(n_cf[, ss == 3, drop = F] * f_cf[, ss == 3, drop = F])
        Nc <- Nc + N_triple
    }
    if(sum(ss == 4) > 0){
        N_quad <- sum(ss == 4) * rowSums(n_cf[, ss == 4, drop = F]) /
            rowSums(n_cf[, ss == 4, drop = F] * f_cf[, ss == 4, drop = F])
        Nc <- Nc + N_quad
    }
    return(Nc)
}


#' Calculate CAI
#'
#' \code{get_cai} calculates Codon Adaptation Index (CAI) of each input CDS
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param rscu rscu table containing CAI weight for each codon. This table could be
#'   generated with \code{est_rscu} or prepared manually.
#' @param level "subfam" (default) or "amino_acid". For which level to determine CAI.
#' @returns a named vector of CAI values
#' @importFrom data.table ':='
#' @importFrom data.table .N
#' @references Sharp PM, Li WH. 1987. The codon Adaptation Index--a measure of directional
#'   synonymous codon usage bias, and its potential applications. Nucleic Acids Res 15:1281-1295.
#' @export
#' @examples
#' # estimate CAI of yeast genes based on RSCU of highly expressed genes
#' heg <- head(yeast_exp[order(-yeast_exp$fpkm), ], n = 500)
#' cf_all <- count_codons(yeast_cds)
#' cf_heg <- cf_all[heg$gene_id, ]
#' rscu_heg <- est_rscu(cf_heg)
#' cai <- get_cai(cf_all, rscu_heg)
#' head(cai)
#' hist(cai)
#'
get_cai <- function(cf, rscu, level = 'subfam'){
    ss <- . <- NULL
    if(!level %in% c('amino_acid', 'subfam')){
      stop('Possible values for `level` are "amino_acid" and "subfam"')
    }
    # exclude single codon sub-family
    rscu <- data.table::as.data.table(rscu)
    rscu[, ss := .N, by = level]
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
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param trna_w tRNA weight for each codon, can be generated with \code{est_trna_weight()}.
#' @returns a named vector of TAI values
#' @references dos Reis M, Savva R, Wernisch L. 2004. Solving the riddle of codon usage
#'   preferences: a test for translational selection. Nucleic Acids Res 32:5036-5044.
#' @export
#' @examples
#' # calculate TAI of yeast genes based on genomic tRNA copy numbers
#' w <- est_trna_weight(yeast_trna_gcn)
#' cf_all <- count_codons(yeast_cds)
#' tai <- get_tai(cf_all, w)
#' head(tai)
#' hist(tai)
#'
get_tai <- function(cf, trna_w){
    # codon frequency per CDS
    cf <- cf[, trna_w$codon, drop = FALSE]

    # tai
    tai <- exp(cf %*% matrix(log(trna_w$w)) / rowSums(cf))
    return(tai[, 1])
}


#' GC contents
#'
#' Calculate GC content of the whole sequences.
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @returns a named vector of GC contents.
#' @export
#' @examples
#' # estimate GC content of yeast genes
#' cf_all <- count_codons(yeast_cds)
#' gc <- get_gc(cf_all)
#' head(gc)
#' hist(gc)
#'
get_gc <- function(cf){
    codon_gc <- sapply(strsplit(colnames(cf), ''), \(x) sum(x %in% c('C', 'G')))
    gc <- cf %*% matrix(codon_gc) / (rowSums(cf) * 3)
    return(gc[, 1])
}


#' GC contents at synonymous 3rd codon positions
#'
#' Calculate GC content at synonymous 3rd codon positions.
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or \code{create_codon_table}.
#' @param level "subfam" (default) or "amino_acid". For which level to determine 
#'   GC content at synonymous 3rd codon positions.
#' @returns a named vector of GC3s values.
#' @importFrom data.table ':='
#' @importFrom data.table .N
#' @references Peden JF. 2000. Analysis of codon usage.
#' @export
#' @examples
#' # estimate GC3s of yeast genes
#' cf_all <- count_codons(yeast_cds)
#' gc3s <- get_gc3s(cf_all)
#' head(gc3s)
#' hist(gc3s)
#'
get_gc3s <- function(cf, codon_table = get_codon_table(), level = 'subfam'){
    aa_code <- ss <- . <- gc3s <- codon <- NULL
    codon_table <- data.table::as.data.table(codon_table)
    codon_table[, ss := .N, by = level]
    codon_table <- codon_table[aa_code != '*' & ss > 1]
    codon_table[, gc3s := substr(codon, 3, 3) %in% c('G', 'C')]

    cf <- cf[, codon_table$codon, drop = FALSE]
    n <- rowSums(cf)
    gc <- rowSums(cf[, codon_table[, codon[gc3s == T]], drop = FALSE])
    return(gc/n)
}


#' GC contents at 4-fold degenerate sites
#'
#' Calculate GC content at synonymous position of codons (using four-fold degenerate sites only).
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @param level "subfam" (default) or "amino_acid". For which level to determine 
#'   GC contents at 4-fold degenerate sites.
#' @returns a named vector of GC4d values.
#' @importFrom data.table ':='
#' @importFrom data.table .N
#' @export
#' @examples
#' # estimate GC4d of yeast genes
#' cf_all <- count_codons(yeast_cds)
#' gc4d <- get_gc4d(cf_all)
#' head(gc4d)
#' hist(gc4d)
#'
get_gc4d <- function(cf, codon_table = get_codon_table(), level = 'subfam'){
    ss <- . <- gc4d <- codon <- NULL
    if(!level %in% c('amino_acid', 'subfam')){
      stop('Possible values for `level` are "amino_acid" and "subfam"')
    }
    codon_table <- data.table::as.data.table(codon_table)
    codon_table[, ss := .N, by = level]
    codon_table <- codon_table[ss == 4]
    codon_table[, gc4d := substr(codon, 3, 3) %in% c('G', 'C')]

    cf <- cf[, codon_table$codon, drop = FALSE]
    n <- rowSums(cf)
    gc <- rowSums(cf[, codon_table[, codon[gc4d == T]], drop = FALSE])
    return(gc/n)
}


#' Fraction of optimal codons (Fop)
#'
#' \code{get_fop} calculates the fraction of optimal codons (Fop) of each CDS.
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param op a character vector of optimal codons. Can be determined automatically by running
#'   \code{est_optimal_codons}.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @param ... other arguments passed to \code{est_optimal_codons}.
#' @returns a named vector of fop values.
#' @references Ikemura T. 1981. Correlation between the abundance of Escherichia coli transfer RNAs
#'   and the occurrence of the respective codons in its protein genes: a proposal for a synonymous
#'   codon choice that is optimal for the E. coli translational system. J Mol Biol 151:389-409.
#' @export
#' @examples
#' # estimate Fop of yeast genes
#' cf_all <- count_codons(yeast_cds)
#' fop <- get_fop(cf_all)
#' head(fop)
#' hist(fop)
#'
get_fop <- function(cf, op = NULL, codon_table = get_codon_table(), ...){
    coef <- qvalue <- codon <- optimal <- NULL
    if(is.null(op)){
        optimal_codons <- est_optimal_codons(cf, codon_table = codon_table, ...)
        op <- optimal_codons[optimal == TRUE, codon]
    }
    rowSums(cf[, op]) / rowSums(cf)
}


#' Mean Codon Stabilization Coefficients
#'
#' \code{get_cscg} calculates Mean Codon Stabilization Coefficients of each CDS.
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param csc table of Codon Stabilization Coefficients as calculated by \code{est_csc()}.
#' @returns a named vector of cscg values.
#' @references Presnyak V, Alhusaini N, Chen YH, Martin S, Morris N, Kline N, Olson S, Weinberg D,
#'   Baker KE, Graveley BR, et al. 2015. Codon optimality is a major determinant of mRNA stability.
#'   Cell 160:1111-1124.
#' @export
#' @examples
#' # estimate CSCg of yeast genes
#' yeast_csc <- est_csc(yeast_cds, yeast_half_life)
#' cf_all <- count_codons(yeast_cds)
#' cscg <- get_cscg(cf_all, csc = yeast_csc)
#' head(cscg)
#' hist(cscg)
#'
get_cscg <- function(cf, csc){
    cf <- cf[, csc$codon, drop = FALSE]
    cp <- cf / rowSums(cf)
    cscg <- cp %*% as.matrix(csc$csc)
    stats::setNames(cscg[, 1], rownames(cscg))
}

#' Deviation from Proportionality
#'
#' \code{get_dp} calculates Deviation from Proportionality of each CDS.
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param host_weights a named vector of tRNA weights for each codon that reflects the relative
#'  availability of tRNAs in the host organism.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @param level "subfam" (default) or "amino_acid". If "subfam", the deviation is calculated at
#'   the codon subfamily level. Otherwise, the deviation is calculated at the amino acid level.
#' @param missing_action Actions to take when no codon of a group were found in a CDS. Options are
#'   "ignore" (default), or "zero" (set codon proportions to 0).
#' @returns a named vector of dp values.
#' @references Chen F, Wu P, Deng S, Zhang H, Hou Y, Hu Z, Zhang J, Chen X, Yang JR. 2020.
#'   Dissimilation of synonymous codon usage bias in virus-host coevolution due to translational
#'   selection. Nat Ecol Evol 4:589-600.
#' @export
#' @examples
#' # estimate DP of yeast genes
#' cf_all <- count_codons(yeast_cds)
#' trna_weight <- est_trna_weight(yeast_trna_gcn)
#' trna_weight <- setNames(trna_weight$w, trna_weight$codon)
#' dp <- get_dp(cf_all, host_weights = trna_weight)
#' head(dp)
#' hist(dp)
#'
get_dp <- function(cf, host_weights, codon_table = get_codon_table(),
                   level = 'subfam', missing_action = 'ignore'){
    aa_code <- NULL # due to NSE notes in R CMD check
    if(!level %in% c('amino_acid', 'subfam')){
      stop('Possible values for `level` are "amino_acid" and "subfam"')
    }
    codon_table <- data.table::as.data.table(codon_table)
    codon_table <- codon_table[!aa_code == '*']
    codon_list <- split(codon_table$codon, codon_table[[level]])
    codon_list <- codon_list[lengths(codon_list) > 1]  # exclude single codon groups

    d <- sapply(codon_list, function(codons){
        # get codon proportions
        cf_grp <- cf[, codons, drop = FALSE]
        codon_prop <- cf_grp / rowSums(cf_grp)
        if(missing_action == 'zero'){
            codon_prop[is.na(codon_prop)] <- 0
        }
        # get relative codon availability in the host organisms
        w <- host_weights[codons]
        rel_avail <- w / sum(w)
        # calculate deviation from proportionality
        sqrt(rowSums(sweep(codon_prop, 2, rel_avail)^2))
    })
    # geometric mean across subfamilies
    # remove NA values due to following causes:
    # 1. no codon of a group were found in a CDS when missing_action = "ignore";
    # 2. exact match to the host tRNA pool. (very rare)
    dp <- exp(rowMeans(log(d), na.rm = TRUE))
}
