#' Calculate ENC
#'
#' \code{get_enc} computes ENC of each CDS
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param codon_table codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @return vector of ENC values, sequence names are used as vector names
#' @export
#' @references
#' * Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87:23-29.
#' * Sun X, Yang Q, Xia X. 2013. An improved implementation of effective number of codons (nc). Mol Biol Evol 30:191-196.
#' @examples
#' # estimate ENC of yeast genes
#' cf_all <- count_codons(yeast_cds)
#' enc <- get_enc(cf_all)
#' head(enc)
#' hist(enc)
get_enc <- function(cf, codon_table = get_codon_table()){
    aa_code <- NULL # due to NSE notes in R CMD check
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
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param rscu rscu table containing CAI weight for each codon. This table could be
#'   generated with `est_rscu` or prepared manually.
#' @returns a named vector of CAI values
#' @importFrom data.table ':='
#' @importFrom data.table .N
#' @references Sharp PM, Li WH. 1987. The codon Adaptation Index--a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Res 15:1281-1295.
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
get_cai <- function(cf, rscu){
    ss <- . <- subfam <- NULL
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
#' @references dos Reis M, Savva R, Wernisch L. 2004. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res 32:5036-5044.
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
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
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
get_gc3s <- function(cf, codon_table = get_codon_table()){
    aa_code <- ss <- . <- subfam <- gc3s <- codon <- NULL
    codon_table[, ss := .N, by = .(subfam)]
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
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
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
get_gc4d <- function(cf, codon_table = get_codon_table()){
    ss <- . <- subfam <- gc4d <- codon <- NULL

    codon_table[, ss := .N, by = .(subfam)]
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
#' @param seqs CDS sequences of all protein-coding genes. One for each gene.
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @returns a named vector of fop values.
#' @references Ikemura T. 1981. Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes: a proposal for a synonymous codon choice that is optimal for the E. coli translational system. J Mol Biol 151:389-409.
#' @export
#' @examples
#' # estimate Fop of yeast genes
#' fop <- get_fop(yeast_cds)
#' head(fop)
#' hist(fop)
#'
get_fop <- function(seqs, codon_table = get_codon_table()){
    coef <- qvalue <- codon <- NULL
    cf <- count_codons(seqs)
    optimal_codons <- est_optimal_codons(seqs, codon_table = codon_table)
    op <- optimal_codons[coef < 0 & qvalue < 0.001, codon]
    rowSums(cf[, op]) / rowSums(cf)
}


#' Mean Codon Stabilization Coefficients
#'
#' \code{get_cscg} calculates Mean Codon Stabilization Coefficients of each CDS.
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param csc table of Codon Stabilization Coefficients as calculated by `est_csc()`.
#' @returns a named vector of cscg values.
#' @references Presnyak V, Alhusaini N, Chen YH, Martin S, Morris N, Kline N, Olson S, Weinberg D, Baker KE, Graveley BR, et al. 2015. Codon optimality is a major determinant of mRNA stability. Cell 160:1111-1124.
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
    cf <- cf[, csc$codon]
    cp <- cf / rowSums(cf)
    cscg <- cp %*% as.matrix(csc$csc)
    stats::setNames(cscg[, 1], rownames(cscg))
}
#' Calculates enc, fop, gc, gc3s, gc4d, cai, tai and cscg of each CDS
#'
#'\code{get_cubar} calculates enc, fop, gc, gc3s, gc4d, cai, tai and cscg of each CDS
#'
get_cubar <- function(seqs, rscu = NULL, trna_w = NULL, csc = NULL){
  seqs_cds_qc <- check_cds(seqs)
  seqs_cf <- count_codons(seqs_cds_qc)
  cf_all <- count_codons(seqs)
  
  enc <- get_enc(seqs_cf)
  fop <- get_fop(seqs)
  gc <- get_gc(cf_all)
  gc3s <- get_gc3s(cf_all)
  gc4d <- get_gc4d(cf_all)

  lengths <- sapply(list(enc, fop, gc, gc3s, gc4d), length)

  # Find the minimum length among columns
  min_length <- min(lengths)

  # Subset each column to the minimum length
  enc <- enc[1:min_length]
  fop <- fop[1:min_length]
  gc <- gc[1:min_length]
  gc3s <- gc3s[1:min_length]
  gc4d <- gc4d[1:min_length]

  # Combine columns
  result <- data.table(enc, fop, gc, gc3s, gc4d)
  
  if (!is.null(rscu)) {
    cai <- get_cai(seqs_cf, rscu)
    cai <- cai[1:min_length]
    result <- cbind(result, cai)
  } 
  
  if (!is.null(trna_w)) {
    tai <- get_tai(seqs_cf, trna_w)
    tai <- tai[1:min_length]
    result <- cbind(result, tai)
  }
  
  if (!is.null(csc)) {
    cscg <- get_cscg(cf_all, csc)
    cscg <- cscg[1:min_length]
    result <- cbind(result, cscg)
  }
  
  result <- as.data.table(result, rn = "seq_id")
  return(result)
}

#' \code{sliding_window_analysis} performs sliding window analysis on CDS and calculates a series of metrics.
#' @param seqs: DNA sequence data represented as a DNAStringSet object, containing DNA sequences of multiple genes.
#' @param window_size: Window size, indicating the length of each sliding window.
#' @param slide_frequency: Slide frequency, indicating the step size between windows.
#' @param csc: Optional parameter to calculate Codon Usage Bias (CUB) related metrics. Default is NULL.
#' @param rscu: Optional parameter to calculate Relative Synonymous Codon Usage (RSCU) related metrics. Default is NULL.
#' @param trna_w: Optional parameter to calculate tRNA adaptation index (tAI) related metrics. Default is NULL.
#' @returns A data frame containing the computed metrics.
sliding_window_analysis <- function(seqs, window_size, slide_frequency, csc = NULL, rscu = NULL, trna_w = NULL) {
  num_genes <- length(seqs)
  result <- DNAStringSet()
  
  for (i in 1:num_genes) {
    gene_data <- as.character(seqs[i])
    num_elements <- nchar(gene_data)
    num_windows <- floor((num_elements - window_size) / slide_frequency) + 1
    
    gene_name <- names(seqs)[i]
    
    if (num_windows <= 0) {
      result <- c(result, DNAStringSet(gene_data))
    } else {
      start_indices <- seq(1, num_elements - window_size + 1, by = slide_frequency)
      end_indices <- start_indices + window_size - 1
      
      windows <- DNAStringSet()
      
      for (j in 1:length(start_indices)) {
        window_seq <- DNAStringSet(substr(gene_data, start_indices[j], end_indices[j]))
        windows <- c(windows, window_seq)
      }
      
      names(windows) <- paste(gene_name, seq_along(windows), sep = "_")
      
      result <- c(result, windows)
    }
  }
  
  result <- DNAStringSet(result)

  result_cf <- count_codons(result)
  gc <- get_gc(result_cf)
  gc3s <- get_gc3s(result_cf)
  gc4d <- get_gc4d(result_cf)
  enc <- get_enc(result_cf)
  
  # To ensure that the fop value is meaningful, when the window size is greater than or equal to 1000, calculate get_fop
  if (window_size >= 1000) {
    fop <- get_fop(result)
  } else {
    fop <- NULL
  }
  
  output <- cbind(gc, gc3s, gc4d, enc, fop)
  
  if (!is.null(csc)) {
    cscg <- get_cscg(result_cf, csc)
    output <- cbind(output, cscg)
  }
  
  if (!is.null(rscu)) {
    cai <- get_cai(result_cf, rscu)
    output <- cbind(output, cai)
  }
  
  if (!is.null(trna_w)) {
    tai <- get_tai(result_cf, trna_w)
    output <- cbind(output, tai)
  }
  
  return(output)
}