#' Calculate effective number of codons (ENC)
#'
#' \code{get_enc} computes the effective number of codons (ENC) for each coding 
#' sequence, which quantifies the degree of codon usage bias. Lower ENC values 
#' indicate stronger bias (fewer codons are used), while higher values indicate 
#' more uniform codon usage.
#'
#' @param cf A matrix of codon frequencies as calculated by \code{count_codons()}.
#'   Rows represent sequences and columns represent codons.
#' @param codon_table A codon table defining the genetic code, derived from 
#'   \code{get_codon_table()} or \code{create_codon_table()}.
#' @param level Character string specifying the analysis level: "subfam" (default, 
#'   analyzes codon subfamilies) or "amino_acid" (analyzes at amino acid level).
#' @return A named numeric vector of ENC values. Names correspond to sequence 
#'   identifiers from the input matrix. ENC values typically range from 20 
#'   (maximum bias) to 61 (uniform usage).
#' @export
#' @references
#' Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87:23-29.
#'
#' Sun X, Yang Q, Xia X. 2013. An improved implementation of effective number of codons (NC). Mol Biol Evol 30:191-196.
#' @examples
#' # Calculate ENC for yeast genes
#' cf_all <- count_codons(yeast_cds)
#' enc <- get_enc(cf_all)
#' head(enc)
#' hist(enc, main = "Distribution of ENC values")
#'
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
        return(rowSums(mx))
    })

    if(nrow(cf) == 1){
        f_cf <- t(f_cf)
        n_cf <- t(n_cf)
    }
    ss <- lengths(codon_list)
    unique_deg <- unique(ss)
    Nc <- sum(ss == 1)
    for (deg in unique_deg[unique_deg > 1]) {
      N_deg <- sum(ss == deg) * rowSums(n_cf[, ss == deg, drop = FALSE]) /
        rowSums(n_cf[, ss == deg, drop = FALSE] * f_cf[, ss == deg, drop = FALSE])
      N_deg[rowSums(n_cf[, ss == deg, drop = FALSE]) == 0] <- sum(ss == deg)* deg
      Nc <- Nc + N_deg
    }
    return(Nc)
}


#' Calculate Codon Adaptation Index (CAI)
#'
#' \code{get_cai} calculates the Codon Adaptation Index (CAI) for each input 
#' coding sequence. CAI measures how similar the codon usage of a gene is to 
#' that of highly expressed genes, serving as an indicator of translational 
#' efficiency. Higher CAI values suggest better adaptation to the translational 
#' machinery.
#'
#' @param cf A matrix of codon frequencies as calculated by \code{count_codons()}.
#'   Rows represent sequences and columns represent codons.
#' @param rscu An RSCU table containing CAI weights for each codon. This table 
#'   should be generated using \code{est_rscu()} based on highly expressed genes, 
#'   or prepared manually with appropriate weight values.
#' @param level Character string specifying the analysis level: "subfam" (default, 
#'   analyzes codon subfamilies) or "amino_acid" (analyzes at amino acid level).
#' @return A named numeric vector of CAI values ranging from 0 to 1. Names 
#'   correspond to sequence identifiers from the input matrix. Values closer 
#'   to 1 indicate higher similarity to highly expressed genes.
#' @importFrom data.table ':='
#' @importFrom data.table .N
#' @references Sharp PM, Li WH. 1987. The codon Adaptation Index--a measure of directional
#'   synonymous codon usage bias, and its potential applications. Nucleic Acids Res 15:1281-1295.
#' @export
#' @examples
#' # Calculate CAI for yeast genes based on RSCU of highly expressed genes
#' heg <- head(yeast_exp[order(-yeast_exp$fpkm), ], n = 500)
#' cf_all <- count_codons(yeast_cds)
#' cf_heg <- cf_all[heg$gene_id, ]
#' rscu_heg <- est_rscu(cf_heg)
#' cai <- get_cai(cf_all, rscu_heg)
#' head(cai)
#' hist(cai, main = "Distribution of CAI values")
#'
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


#' Calculate tRNA Adaptation Index (TAI)
#'
#' \code{get_tai} calculates the tRNA Adaptation Index (TAI) for each coding 
#' sequence, which measures how well codon usage matches tRNA availability in 
#' the cell. Higher TAI values indicate better adaptation to the tRNA pool, 
#' suggesting more efficient translation.
#'
#' @param cf A matrix of codon frequencies as calculated by \code{count_codons()}. 
#'   Note: Start codons should be removed from sequences before analysis to 
#'   avoid bias from universal start codon usage.
#' @param trna_w A table of tRNA weights for each codon, generated using 
#'   \code{est_trna_weight()}. These weights reflect relative tRNA availability.
#' @param w_format Character string specifying the format of tRNA weights: 
#'   "cubar" (default, weights from cubar package) or "tAI" (weights from 
#'   the tAI package format).
#' @return A named numeric vector of TAI values. Names correspond to sequence 
#'   identifiers from the input matrix. Values range from 0 to 1, with higher 
#'   values indicating better adaptation to tRNA availability.
#' @references dos Reis M, Savva R, Wernisch L. 2004. Solving the riddle of codon usage
#'   preferences: a test for translational selection. Nucleic Acids Res 32:5036-5044.
#' @export
#' @examples
#' # calculate TAI of yeast genes based on genomic tRNA copy numbers
#' w <- est_trna_weight(yeast_trna_gcn)
#' yeast_cds_qc <- check_cds(yeast_cds)
#' cf <- count_codons(yeast_cds_qc)
#' tai <- get_tai(cf, w)
#' head(tai)
#' hist(tai)
#'
get_tai <- function(cf, trna_w, w_format = "cubar"){
    
    if(!w_format %in% c('cubar', 'tAI')){
      stop('Possible values for `w_format` are "cubar" and "tAI"')
    }
    if(w_format == "cubar"){
      # codon frequency per CDS
      cf <- cf[, trna_w$codon, drop = FALSE]
      # tai
      tai <- exp(cf %*% matrix(log(trna_w$w)) / rowSums(cf))
    }else if(w_format == "tAI"){
      codon_table <- get_codon_table()
      valid_codons <- setdiff(codon_table$codon, c('TAA', 'TAG', 'TGA', 'ATG'))
      tai <- exp(cf[, valid_codons] %*% matrix(log(trna_w)) / rowSums(cf[, valid_codons]))
    }
    return(tai[, 1])
}


#' Calculate GC content of coding sequences
#'
#' \code{get_gc} calculates the overall GC content (percentage of guanine and 
#' cytosine nucleotides) for each coding sequence.
#'
#' @param cf A matrix of codon frequencies as calculated by \code{count_codons()}.
#'   Rows represent sequences and columns represent codons.
#' @return A named numeric vector of GC content values (ranging from 0 to 1). 
#'   Names correspond to sequence identifiers from the input matrix.
#' @export
#' @examples
#' # Calculate GC content for yeast genes
#' cf_all <- count_codons(yeast_cds)
#' gc <- get_gc(cf_all)
#' head(gc)
#' hist(gc, main = "Distribution of GC content")
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
#' @returns a named vector of GC3s values. The names of the elements correspond to the sequence names.
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
#' @returns a named vector of GC4d values. The names of the elements correspond to the sequence names.
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


#' Calculate fraction of optimal codons (Fop)
#'
#' \code{get_fop} calculates the fraction of optimal codons (Fop) for each 
#' coding sequence, which represents the proportion of codons that are 
#' considered optimal for translation efficiency. Higher Fop values suggest 
#' stronger selection for optimal codon usage.
#'
#' @param cf A matrix of codon frequencies as calculated by \code{count_codons()}.
#'   Rows represent sequences and columns represent codons.
#' @param op A character vector specifying which codons are considered optimal. 
#'   If not provided, optimal codons will be determined automatically using 
#'   \code{est_optimal_codons()}.
#' @param codon_table A codon table defining the genetic code, derived from 
#'   \code{get_codon_table()} or \code{create_codon_table()}.
#' @param ... Additional arguments passed to \code{est_optimal_codons()} when 
#'   optimal codons are determined automatically.
#' @return A named numeric vector of Fop values (ranging from 0 to 1). Names 
#'   correspond to sequence identifiers from the input matrix. Higher values 
#'   indicate greater usage of optimal codons.
#' @references Ikemura T. 1981. Correlation between the abundance of Escherichia coli transfer RNAs
#'   and the occurrence of the respective codons in its protein genes: a proposal for a synonymous
#'   codon choice that is optimal for the E. coli translational system. J Mol Biol 151:389-409.
#' @export
#' @examples
#' # Calculate Fop for yeast genes (optimal codons determined automatically)
#' cf_all <- count_codons(yeast_cds)
#' fop <- get_fop(cf_all)
#' head(fop)
#' hist(fop, main = "Distribution of Fop values")
#'
#'
get_fop <- function(cf, op = NULL, codon_table = get_codon_table(), ...){
  coef <- qvalue <- codon <- optimal <- amino_acid <- N <- . <- NULL
  excluded_codon <- codon_table[, .(codon = codon, .N), by = .(amino_acid)][(N == 1 | amino_acid == "*")]$codon
  cf <- cf[, !colnames(cf) %in% excluded_codon]
  codon_table <- codon_table[!codon %in% excluded_codon,]
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
#' @returns a named vector of cscg values. The names of the elements correspond to the sequence names.
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
#' @returns a named vector of dp values. The names of the elements correspond to the sequence names.
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

#' Amino Acid Usage
#'
#' Calculate Amino Acid Usage Frequencies of each CDS.
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @returns a matrix of amino acid frequencies for each CDS.
#' Each row corresponds to a sequence, and each column represents an amino acid.
#' @export
#' @examples
#' # estimate amino acid frequencies of yeast CDSs
#' cf_all <- count_codons(yeast_cds)
#' aau_gene <- get_aau(cf_all)
#' head(aau_gene)
#' 

get_aau <- function(cf, codon_table = get_codon_table()){
  aa_code <- id <- codon <- count <- count_codon <- aa_code <-  amino_acid <- . <- NULL # due to NSE notes in R CMD check
  codon_table <- data.table::as.data.table(codon_table)
  codon_table <- codon_table[amino_acid != '*']
  aa_mat <- sapply(split(codon_table$codon, codon_table$amino_acid), \(x) rowSums(cf[, x, drop = FALSE]))
  if(!is.matrix(aa_mat)){
    aa_mat <- t(as.matrix(aa_mat))
  }
  aau_mat <- aa_mat/rowSums(aa_mat)
  return(aau_mat)
}
