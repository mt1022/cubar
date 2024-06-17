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
get_fop <- function(seqs_cds_qc, codon_table = get_codon_table()){
    coef <- qvalue <- codon <- NULL
    cf <- count_codons(seqs_cds_qc)
    optimal_codons <- est_optimal_codons(seqs_cds_qc, codon_table = codon_table)
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
#' \code{get_cubar} calculates enc, fop, gc, gc3s, gc4d, cai, tai and cscg of each CDS
#'
#' @param rscu table containing CAI weight for each codon. This table could be
#'   generated with `est_rscu` or prepared manually.
#' @param trna_w tRNA weight for each codon, can be generated with `est_trna_weight()`.
#' @param csc table of Codon Stabilization Coefficients as calculated by `est_csc()`.
#'
get_cubar <- function(seqs, rscu = NULL, trna_w = NULL, csc = NULL){
  seqs_cds_qc <- check_cds(seqs)
  seqs_cf <- count_codons(seqs_cds_qc)

  enc <- get_enc(seqs_cf)
  fop <- get_fop(seqs_cds_qc)
  gc <- get_gc(seqs_cf)
  gc3s <- get_gc3s(seqs_cf)
  gc4d <- get_gc4d(seqs_cf)

  lengths <- sapply(list(enc, fop, gc, gc3s, gc4d), length)

  # Combine columns
  result <- cbind(enc, fop, gc, gc3s, gc4d)

  if (!is.null(rscu)) {
    cai <- get_cai(seqs_cf, rscu)
    result <- cbind(result, cai)
  }

  if (!is.null(trna_w)) {
    tai <- get_tai(seqs_cf, trna_w)
    result <- cbind(result, tai)
  }

  if (!is.null(csc)) {
    cscg <- get_cscg(seqs_cf, csc)
    result <- cbind(result, cscg)
  }

  return(result)
}


#' Plots the results obtained by get_cubar()
#'
#' \code{plot_cubar} displays the parameters calculated by get_cubar() in a histogram
#' and visualizes the distribution of all numerical variables for each parameter
#'
#' @param df the parameters calculated by get_cubar().
#' @param rscu table containing CAI weight for each codon. This table could be
#'   generated with `est_rscu` or prepared manually.
#' @param trna_w tRNA weight for each codon, can be generated with `est_trna_weight()`.
#' @param csc table of Codon Stabilization Coefficients as calculated by `est_csc()`.
#' @returns A data frame containing the computed metrics.
#' @examples
#' df_1 <- get_cubar(yeast_cds)
#' plot_cubar(df_1)
#'
#' df_2 <- get_cubar(yeast_cds,rscu = rscu_heg,trna_w = trna_w, csc = yeast_csc)
#' plot_cubar(df_2)

library(ggplot2)
library(tidyr)
plot_cubar <- function(df){
    # Convert data frame from wide format to long format
    long_df <- pivot_longer(df, cols = everything(), names_to = "variable", values_to = "value")

    # Use ggplot2 to plot the histogram and facet_wrap to create a facet for each variable
    p <- ggplot(long_df, aes(x = value)) +
        geom_histogram(bins = 30) + #Adjust the number of bins as needed
        facet_wrap(~ variable, scales = "free") +
        labs(x = "values", y = "Number of genes") +
        theme_minimal()

    print(p)
}


#' Analyzes gene sequences with a sliding window to get parameters for each segment.
#'
#' \code{get_windowAnalysis} performs sliding window analysis on CDS
#' and calculates enc, fop, gc, gc3s, gc4d, cai, tai and cscg of each CDS.
#'
#' @param seqs DNA sequence data represented as a DNAStringSet object, containing DNA sequences of multiple genes.
#' @param window_size indicates the length of each sliding window.
#' @param slide_frequency indicates the step size between windows.
#' @param rscu table containing CAI weight for each codon. This table could be
#'   generated with `est_rscu` or prepared manually.
#' @param trna_w tRNA weight for each codon, can be generated with `est_trna_weight()`.
#' @param csc table of Codon Stabilization Coefficients as calculated by `est_csc()`.
#' @returns A data frame containing the computed metrics.

get_windowAnalysis <- function(seqs, window_size, slide_frequency, csc = NULL, rscu = NULL, trna_w = NULL) {
  seqs <- check_cds(seqs)
  num_genes <- length(seqs)
  result <- vector("list", num_genes)

  for (i in 1:num_genes) {
    gene_data <- as.character(seqs[i])
    num_elements <- nchar(gene_data)
    num_windows <- floor((num_elements - window_size) / slide_frequency) + 1

    if (num_windows <= 0) {
      result[[i]] <- DNAStringSet(gene_data)
    } else {
      start_indices <- seq(1, num_elements - window_size + 1, by = slide_frequency)
      end_indices <- start_indices + window_size - 1

      windows <- vector("list", length(start_indices))
      for (j in seq_along(start_indices)) {
        window_seq <- substr(gene_data, start_indices[j], end_indices[j])
        windows[[j]] <- DNAString(window_seq)  # Ensure each element is a DNAString object
      }
      names(windows) <- paste(names(seqs)[i], seq_along(windows), sep = "_")
      result[[i]] <- DNAStringSet(do.call(c, windows))  # Combine using DNAStringSet constructor
    }
  }

    result <- do.call(c, result)
  result <- DNAStringSet(result)

  result_cf <- count_codons(result)
  gc <- get_gc(result_cf)
  gc3s <- get_gc3s(result_cf)
  gc4d <- get_gc4d(result_cf)
  enc <- get_enc(result_cf)

  output <- data.frame(gc, gc3s, gc4d, enc)

  # Calculate FOP only if the number of genes is greater than or equal to 1000
  if (num_genes >= 1000) {
    fop <- get_fop(result)
    output <- cbind(output, fop)
  }

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

#' Plots the results obtained by get_windowAnalysis()
#'
#' \code{plot_window} displays the parameters calculated by get_windowAnalysis() in a line chart.
#'
#' @param df the parameters calculated by get_windowAnalysis().
#' @returns a line chart containing the computed metrics.
#' @examples
#' df_1 <- get_windowAnalysis(yeast_cds[1:10],30,30,rscu = rscu_heg,trna_w = trna_w, csc = yeast_csc )
#' plot_window(df_1)
#'
#' df_2 <- get_windowAnalysis(yeast_cds,300,30 )
#' plot_window(df_2)
#'
#' df_3 <- get_windowAnalysis(yeast_cds[1:20],300,30)
#' plot_window(df_3)

library(ggplot2)
library(dplyr)
library(tidyr)
plot_window <- function(df) {

  df$fragment <- rownames(df)
  df$fragment_number <- NA

 # Find the fragments that end with a number and update fragment_number
  numeric_fragments <- grepl(".*[_-]\\d+$", df$fragment)
  df$fragment_number[numeric_fragments] <- as.numeric(sub(".*[_-](\\d+)$", "\\1", df$fragment[numeric_fragments]))

  # For fragments that do not end with a number, append '_1' and update fragment_number
  non_numeric_fragments <- !numeric_fragments
  df$fragment[non_numeric_fragments] <- paste0(df$fragment[non_numeric_fragments], "_1")
  df$fragment_number[non_numeric_fragments] <- 1

  df <- df[!is.na(df$fragment_number), ]

  # If no rows remain after filtering, stop the function and return a message
  if (nrow(df) == 0) {
    stop("All fragment identifiers are non-standard and have been removed.")
  }

  df <- df %>% select(-fragment)

  df_long <- pivot_longer(df, cols = -fragment_number, names_to = "gene", values_to = "value")

  if (!is.numeric(df_long$value)) {
    stop("The 'value' column must be numeric.")
  }

  df_means <- df_long %>%
    group_by(fragment_number, gene) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = 'drop')

  p <- ggplot(df_means, aes(x = fragment_number, y = value)) +
    geom_line() +
    facet_wrap(~ gene, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Fragment Number", y = "Average Value")
  return(p)

}

