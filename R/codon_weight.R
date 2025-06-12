#' Estimate RSCU
#'
#' \code{est_rscu} returns the RSCU value of codons
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param weight a vector of the same length as \code{seqs} that gives different weights
#'   to CDSs when count codons. for example, it could be gene expression levels.
#' @param pseudo_cnt pseudo count to avoid dividing by zero. This may occur when
#'   only a few sequences are available for RSCU calculation.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @param level "subfam" (default) or "amino_acid". For which level to determine RSCU.
#' @returns a data.table of codon info. RSCU values are reported in the last column.
#' @importFrom data.table ':='
#' @references Sharp PM, Tuohy TM, Mosurski KR. 1986. Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucleic Acids Res 14:5125-5143.
#' @export
#'
#' @examples
#' # compute RSCU of all yeast genes
#' cf_all <- count_codons(yeast_cds)
#' est_rscu(cf_all)
#'
#' # compute RSCU of highly expressed (top 500) yeast genes
#' heg <- head(yeast_exp[order(-yeast_exp$fpkm), ], n = 500)
#' cf_heg <- count_codons(yeast_cds[heg$gene_id])
#' est_rscu(cf_heg)
#'
est_rscu <- function(cf, weight = 1, pseudo_cnt = 1, codon_table = get_codon_table(),
                     level = 'subfam'){
    aa_code <- cts <- codon <- . <- rscu <- prop <- NULL # due to NSE notes in R CMD check
    if(!level %in% c('amino_acid', 'subfam')){
      stop('Possible values for `level` are "amino_acid" and "subfam"')
    }
    codon_freq <- colSums(cf * weight)
    codon_table <- codon_table[aa_code != '*']
    codon_table[, cts := codon_freq[codon]]
    codon_table[, `:=`(
        prop = (cts + pseudo_cnt) / sum(cts + pseudo_cnt),
        w_cai = (cts + pseudo_cnt) / max(cts + pseudo_cnt)), by = level]
    codon_table[, rscu := prop / mean(prop), by = level]
    return(codon_table[])
}


#' Generate codon-anticodon pairing relationship
#'
#' \code{ca_pairs} show possible codon-anticodons pairings
#'
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @param domain "Eukarya" (default), "Bacteria" or "Archaea". Determine based on the type of organism.
#' @param plot FALSE (default) or TRUE. Whether to keep the columns required for plotting.
#' @returns a data.table of codon-anticodon pairing information
#' @importFrom data.table ':='
#' @importFrom rlang .data
#' @export
#' @examples
#' # get possible codon and anticodon pairings for the vertebrate mitochondrial genetic code
#' ctab <- get_codon_table(gcid = '2')
#' pairing <- ca_pairs(ctab)
#' head(pairing)
#'
ca_pairs <- function(codon_table = get_codon_table(), domain = "Eukarya", plot = FALSE){
    anticodon <- codon <- codon_b1 <- codon_b2 <- codon_b3 <- amino_acid <- NULL # due to NSE notes in R CMD check
    . <- aa_code <- base_codon <- base_anti <- type <- NULL
    anticodon_aa <- codon_aa <- i.aa_code <- i.amino_acid <- NULL
    
    if(!domain %in% c("Eukarya", "Bacteria", "Archaea")){
      stop("Unknown domain: ", domain)
    }
    codon_table <- data.table::copy(codon_table)
    
    if(nrow(codon_table) < 64){
      message("Warning: The input codon table is incomplete. The missing codons are filled in as termination codons.")
      codon_table <- codon_table[names(Biostrings::GENETIC_CODE), on = .(codon)]
      codon_table[is.na(aa_code), `:=`(
        aa_code = '*', amino_acid = '*', subfam = paste('*', substr(codon, 1, 2), sep = '_'))]
    }
    codon_table[, anticodon := as.character(rev_comp(codon_table$codon))]
    codon_table[, c('codon_b1', 'codon_b2', 'codon_b3') := data.table::tstrsplit(codon, '')]
    bases <- c('T', 'C', 'A', 'G')
    codon_table[, codon_b1 := factor(codon_b1, levels = bases)]
    codon_table[, `:=`(
        codon_b1 = factor(codon_b1, levels = bases),
        codon_b2 = factor(codon_b2, levels = bases),
        codon_b3 = factor(codon_b3, levels = rev(bases))
    )]
    codon_table <- codon_table[order(codon_b1, codon_b2, codon_b3)]
    ca_pairs <- codon_table[, {
        wc <- data.table::data.table(
            type = 'WC',
            base_codon = codon_b3[1:3],
            base_anti = codon_b3[c(1:3)],
            codon = codon[1:3], anticodon = anticodon[1:3])
        # type: anticodon base + corresponding codon base
        # note: base_anti here is used for visualization, not referring to 1st base of anticodon
        ih <- data.table::data.table(
            type = c('IU', 'IC', 'IA'),
            base_codon = codon_b3[c(4, 3, 2)],
            base_anti = codon_b3[c(4, 4, 4)],
            codon = codon[c(4, 3, 2)], anticodon = anticodon[4])
        gu <- data.table::data.table(
            type = c('GU', 'UG'),
            base_codon = codon_b3[c(4, 1)],
            base_anti = codon_b3[c(3, 2)],
            codon = codon[c(4, 1)], anticodon = anticodon[c(3, 2)])
        rbind(wc, ih, gu)
    }, by = .(codon_b1, codon_b2)]
    stop_codons <- codon_table[aa_code == '*', codon]
    # no pairing to stop codons or anitcodon corresponding to stop codons
    ca_pairs <- ca_pairs[!codon %in% stop_codons]
    ca_pairs <- ca_pairs[!anticodon %in% as.character(rev_comp(stop_codons))]
    # only wobble among synonymous codons and anticodons
    ca_pairs[codon_table, codon_aa := i.aa_code, on = .(codon)]
    ca_pairs[codon_table, amino_acid := i.amino_acid, on = .(codon)]
    ca_pairs[codon_table, anticodon_aa := i.aa_code, on = .(anticodon)]
    ca_pairs <- ca_pairs[codon_aa == anticodon_aa]
    
    if(domain == "Bacteria"){
      ca_pairs <- rbind(ca_pairs, data.table::data.table(codon_b1 = "A", codon_b2 = "T", type = "LA", base_codon = "A", base_anti = "G",
                                                         codon = "ATA", anticodon = "CAT", codon_aa = "I", amino_acid = "Ile", anticodon_aa = "I"))
    }else if(domain == "Archaea"){
      ca_pairs <- rbind(ca_pairs, data.table::data.table(codon_b1 = "A", codon_b2 = "T", type = "agmA", base_codon = "A", base_anti = "G",
                                                         codon = "ATA", anticodon = "CAT", codon_aa = "I", amino_acid = "Ile", anticodon_aa = "I"))
    }
    ca_pairs <- ca_pairs[order(factor(ca_pairs$codon, levels = codon_table$codon))]
    if(plot){
      return(ca_pairs)
    }else{
      return(ca_pairs[, .(type, codon, anticodon, amino_acid)])
    }
    
}


#' Plot codon-anticodon pairing relationship
#'
#' \code{plot_ca_pairs} show possible codon-anticodons pairings
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @param ca_pairs a table of codon-anticodon pairing derived from \code{ca_pairs}
#' @returns a plot on possible codon-anticodons pairings
#' @importFrom data.table ':='
#' @importFrom rlang .data
#' @export
#' @examples
#' # plot possible codon and anticodon pairings for the vertebrate mitochondrial genetic code
#' ctab <- get_codon_table(gcid = '2')
#' pairing <- ca_pairs(ctab, plot = TRUE)
#' plot_ca_pairs(ctab, pairing)
#'
#' # plot possible codon and anticodon pairings for the standard genetic code in bacteria
#' plot_ca_pairs(ca_pairs = ca_pairs(domain = "Bacteria", plot = TRUE))
#' 
plot_ca_pairs <- function(codon_table = get_codon_table(), ca_pairs = ca_pairs){
  anticodon <- codon <- codon_b1 <- codon_b2 <- codon_b3 <- . <- NULL # due to NSE notes in R CMD check
  aa_code <- base_codon <- base_anti <- type <- NULL
  codon_table <- data.table::copy(codon_table)
  
  if(nrow(codon_table) < 64){
    message("Warning: The input codon table is incomplete. The missing codons are filled in as termination codons.")
    codon_table <- codon_table[names(Biostrings::GENETIC_CODE), on = .(codon)]
    codon_table[is.na(aa_code), `:=`(aa_code = '*', amino_acid = '*', subfam = paste('*', substr(codon, 1, 2), sep = '_'))]
  }
  codon_table[, anticodon := as.character(rev_comp(codon_table$codon))]
  codon_table[, c('codon_b1', 'codon_b2', 'codon_b3') := data.table::tstrsplit(codon, '')]
  bases <- c('T', 'C', 'A', 'G')
  codon_table[, codon_b1 := factor(codon_b1, levels = bases)]
  codon_table[, `:=`(
    codon_b1 = factor(codon_b1, levels = bases),
    codon_b2 = factor(codon_b2, levels = bases),
    codon_b3 = factor(codon_b3, levels = rev(bases))
  )]
  codon_table <- codon_table[order(codon_b1, codon_b2, codon_b3)]
  p <- ggplot2::ggplot(codon_table, ggplot2::aes(y = .data$codon_b3)) +
    ggplot2::geom_segment(
      data = ca_pairs,
      mapping = ggplot2::aes(
        x = 1, xend = 2,
        y = base_codon, yend = base_anti,
        color = type)) +
    ggplot2::geom_label(ggplot2::aes(x = 1, label = .data$codon)) +
    ggplot2::geom_label(ggplot2::aes(x = 2, label = .data$anticodon)) +
    ggplot2::geom_label(ggplot2::aes(x = 0.4, label = .data$aa_code)) +
    ggplot2::scale_color_brewer(palette = 'Dark2') +
    ggplot2::facet_grid(codon_b2 ~ codon_b1, scales = 'free_y') +
    ggplot2::scale_x_continuous(
      limits = c(0.15, 2.5), breaks = c(0.4, 1, 2),
      labels = c('AA', 'Codon', 'Anti')) +
    ggplot2::labs(x = NULL, y = NULL, color = NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   strip.text = ggplot2::element_blank())
  print(p)
}

#' get tRNA gene copy number from GtRNADB
#'
#' \code{trna_gcn} get tRNA gene copy number from GtRNADB
#' @param trna_seq a fasta file of tRNA sequences from GtRNADB
#' @returns a table of tRNA gene copy number for each anticodon
#' @export
#' @examples
#' # get tRNA gene copy number for yeast
#' trna_gcn <- trna_gcn(yeast_trna)
#' 
trna_gcn <- function(trna_seq){
  trna_copy <- sub('.*tRNA-(.*?)-\\d.*', '\\1', names(trna_seq))
  trna_copy <- trna_copy[!trna_copy %in% c("fMet-CAT", "iMet-CAT", "Und-NNN")]
  trna_gcn <- table(as.vector(trna_copy))
  return(trna_gcn)
}


#' Estimate tRNA weight w
#'
#' \code{est_trna_weight} compute the tRNA weight per codon for TAI calculation.
#' This weight reflects relative tRNA availability for each codon.
#'
#' @param trna_level, named vector of tRNA level (or gene copy numbers), one value for each anticodon.
#'   vector names are anticodons.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or \code{create_codon_table}.
#' @param domain "Eukarya" (default), "Bacteria" or "Archaea". Determine based on the type of organism. 
#' The inferred wobble values strength are different for each domain of life. Specify either the parameter "domain" or "s".
#' @param s list of non-Waston-Crick pairing panelty. Specify either the parameter "domain" or "s".
#' @returns data.table of tRNA expression information.
#' @importFrom data.table ':='
#' @references dos Reis M, Savva R, Wernisch L. 2004. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res 32:5036-5044.
#' @export
#' @examples
#' # estimate codon tRNA weight for yeast
#' yeast_trna_w <- est_trna_weight(yeast_trna_gcn)
#' print(yeast_trna_w)
#' 
est_trna_weight <- function(trna_level, codon_table = get_codon_table(), domain = "Eukarya",
                            s = NULL){
    anticodon <- aa_code <- ac_level <- penality <- amino_acid <- NULL # due to NSE notes in R CMD check
    i.values <- . <- ind <- codon <- trna_id <- type <- W <- i.W <- w <- NULL # due to NSE notes in R CMD check
    
    codon_table1 <- data.table::copy(codon_table)
    codon_table[, anticodon := as.character(Biostrings::reverseComplement(
      Biostrings::DNAStringSet(codon_table$codon)))]
    codon_table <- codon_table[aa_code != '*']
    if(domain %in% c("Bacteria", "Archaea")){
      codon_table <- rbind(codon_table, data.table::data.table(aa_code = "I", amino_acid = "Ile", codon = "ATA",
                                                               subfam = "Ile_AT", anticodon = "CAT"))
    }
    codon_table[, trna_id := paste(amino_acid, anticodon, sep = "-")]
    codon_table[, ac_level := trna_level[trna_id]]
    codon_table[is.na(ac_level), ac_level := 0]
    
    # Use the original codon table as input to avoid warnings due to the removal of stop codons
    ca_pairs <- ca_pairs(codon_table = codon_table1, domain = domain)
    
    if(is.null(s)){
      s <- switch(domain,
                  Eukarya  = list(WC = 0, IU = 0, IC = 0.4659, IA = 0.9075, GU = 0.7861, UG = 0.6295),
                  Bacteria = list(WC = 0, IU = 0, IC = 0.4211, IA = 0.8773, GU = 0.6294, UG = 0.698, LA = 0.7309),
                  Archaea  = list(WC = 0, IU = 0, IC = 0.3774, IA = 0.5015, GU = 0.3898, UG = 0.4363, agmA = 0.6453),
                  stop("Unknown domain: ", domain)
      )
    }else{
      if(!"WC" %in% names(s)){
        s$WC <- 0
      }
      missing_types <- setdiff(unique(ca_pairs$type), names(s))
      if(length(missing_types) > 0){
        stop("Missing selection penalty for codon pairs: ", paste(missing_types, collapse = " "))
      }
    }

    s <- utils::stack(s)
    ca_pairs[s, penality := i.values, on = .(type = ind)]
    ca_pairs[, trna_id := paste(amino_acid, anticodon, sep = "-")]
    ca_pairs <- ca_pairs[trna_id %in% names(trna_level)]
    ca_pairs[, ac_level := trna_level[trna_id]]

    dtt_W <- ca_pairs[, .(W = sum(ac_level * (1 - penality))), by = .(codon)]
    codon_table[dtt_W, W := i.W, on = .(codon)]
    codon_table[, w := W/max(W, na.rm = TRUE)]

    # using geometric mean of w for rare cases that a codon has no matching tRNA
    # probably due to incomplete tRNA annotation
    w0 <- codon_table[!is.na(w), w]
    mean_w <- exp(mean(log(w0)))
    codon_table[is.na(w), w := mean_w]
    return(codon_table)
}


#' Estimate optimal codons
#'
#' \code{est_optimal_codons} determine optimal codon of each codon family with binomial regression.
#'   Usage of optimal codons should correlate negatively with enc.
#'
#' @param cf matrix of codon frequencies as calculated by \code{count_codons()}.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or \code{create_codon_table}.
#' @param level "subfam" (default) or "amino_acid". For which level to determine optimal codons.
#' @param gene_score a numeric vector of scores for genes. The order of values should match with
#'   gene orders in the codon frequency matrix. The length of the vector should be equal to the
#'   number of rows in the matrix. The scores could be gene expression levels (RPKM or TPM) that are
#'   optionally log-transformed (for example, with \code{log1p}). The opposite of ENC will be used by
#'   default if \code{gene_score} is not provided.
#' @param fdr false discovery rate used to determine optimal codons.
#' @returns data.table of optimal codons.
#' @importFrom data.table ':='
#' @export
#' @examples
#' # perform binomial regression for optimal codon estimation
#' cf_all <- count_codons(yeast_cds)
#' codons_opt <- est_optimal_codons(cf_all)
#' codons_opt <- codons_opt[optimal == TRUE]
#' codons_opt
#'
est_optimal_codons <- function(cf, codon_table = get_codon_table(), level = 'subfam',
                               gene_score = NULL, fdr = 0.001){
    aa_code <- . <- codon <- NULL # due to NSE notes in R CMD check
    qvalue <- pvalue <- optimal <- coef <- NULL
    if(!level %in% c('amino_acid', 'subfam')){
        stop('Possible values for `level` are "amino_acid" and "subfam"')
    }
    if(is.null(gene_score)){
        enc <- get_enc(cf, codon_table = codon_table)
        # larger values correlate with stronger bias like gene expression levels
        gene_score <- -enc
    }

    # exclude stop codons
    codon_table <- codon_table[aa_code != '*']
    cf <- cf[, colnames(cf) %in% codon_table$codon]

    # regression analysis for each codon sub-family
    binreg <- lapply(split(codon_table$codon, f = codon_table[[level]]), function(x){
        cf_grp <- cf[, x, drop = FALSE]
        if(ncol(cf_grp) == 1){
            data.table::data.table(
                codon = colnames(cf_grp), coef = 0, se = 0, zvalue = 0, pvalue = 0)
        }else{
            total <- rowSums(cf_grp)
            res <- apply(cf_grp, 2, function(x){
                x
                fit <- stats::glm(cbind(x, total - x) ~ gene_score, family = 'binomial')
                summary(fit)$coefficients[-1, ]
            })
            res <- data.table::as.data.table(
                as.data.frame(t(res)), keep.rownames = 'codon')
            data.table::setnames(res, c('codon', 'coef', 'se', 'zvalue', 'pvalue'))
        }
    })
    bingreg <- data.table::rbindlist(binreg, idcol = level)
    bingreg[, qvalue := stats::p.adjust(pvalue, method = 'BH')]
    bingreg <- codon_table[bingreg, on = c('codon', level)]
    bingreg[, c('se', 'zvalue') := NULL]
    bingreg[, optimal := coef > 0 & qvalue < fdr]
    return(bingreg)
}


#' Estimate Codon Stabilization Coefficient
#'
#' \code{get_csc} calculate codon occurrence to mRNA stability correlation coefficients (Default to Pearson's).
#'
#' @param seqs CDS sequences of all protein-coding genes. One for each gene.
#' @param half_life data.frame of mRNA half life (gene_id & half_life are column names).
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or \code{create_codon_table}.
#' @param cor_method method name passed to `cor.test` used for calculating correlation coefficients.
#' @returns data.table of optimal codons.
#' @references Presnyak V, Alhusaini N, Chen YH, Martin S, Morris N, Kline N, Olson S, Weinberg D, Baker KE, Graveley BR, et al. 2015. Codon optimality is a major determinant of mRNA stability. Cell 160:1111-1124.
#' @export
#' @examples
#' # estimate yeast mRNA CSC
#' est_csc(yeast_cds, yeast_half_life)
#'
est_csc <- function(seqs, half_life, codon_table = get_codon_table(), cor_method = 'pearson'){
    aa_code <- codon <- NULL # due to NSE notes in R CMD check
    non_stop_codons <- codon_table[aa_code != '*', codon]
    cf <- count_codons(seqs)

    common_genes <- intersect(half_life$gene_id, rownames(cf))
    half_life <- half_life[half_life$gene_id %in% common_genes, ]
    cf <- cf[half_life$gene_id, non_stop_codons]
    cp <- cf / rowSums(cf)  # codon proportions

    csc <- apply(cp, 2, function(prop){
        res <- stats::cor.test(prop, half_life$half_life, method = cor_method)
        c(corr = res$estimate, pvalue = res$p.value)
    })
    csc <- data.table::data.table(t(csc), keep.rownames = TRUE)
    data.table::setnames(csc, c('codon', 'csc', 'pvalue'))
    return(csc)
}