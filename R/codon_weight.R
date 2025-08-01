#' Estimate Relative Synonymous Codon Usage (RSCU)
#'
#' \code{est_rscu} calculates the Relative Synonymous Codon Usage (RSCU) values 
#' for codons, which quantify the bias in synonymous codon usage. RSCU values 
#' indicate whether a codon is used more (>1) or less (<1) frequently than 
#' expected under uniform usage within its synonymous group.
#'
#' @param cf A matrix of codon frequencies as calculated by \code{count_codons()}.
#'   Rows represent sequences and columns represent codons.
#' @param weight A numeric vector of the same length as the number of sequences 
#'   in \code{cf}, providing different weights for sequences when calculating 
#'   codon frequencies. For example, gene expression levels. Default is 1 
#'   (equal weights).
#' @param pseudo_cnt Numeric pseudo count added to avoid division by zero when 
#'   few sequences are available for RSCU calculation (default: 1).
#' @param codon_table A codon table defining the genetic code, derived from 
#'   \code{get_codon_table()} or \code{create_codon_table()}.
#' @param level Character string specifying the analysis level: "subfam" (default, 
#'   analyzes codon subfamilies) or "amino_acid" (analyzes at amino acid level).
#' @param incl_stop Logical. Whether to include RSCU values for stop codons 
#'   in the output (default: FALSE).
#' @return A data.table containing the codon table with additional columns for 
#'   RSCU analysis: usage frequency counts (cts), frequency proportions (prop), 
#'   CAI weights (w_cai), and RSCU values (rscu). The table includes amino acid 
#'   codes, full amino acid names, codons, and subfamily classifications.
#' @importFrom data.table ':='
#' @references Sharp PM, Tuohy TM, Mosurski KR. 1986. Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucleic Acids Res 14:5125-5143.
#' @export
#' @examples
#' # Calculate RSCU for all yeast genes
#' cf_all <- count_codons(yeast_cds)
#' rscu_all <- est_rscu(cf_all)
#' head(rscu_all)
#'
#' # Calculate RSCU for highly expressed genes (top 500)
#' heg <- head(yeast_exp[order(-yeast_exp$fpkm), ], n = 500)
#' cf_heg <- count_codons(yeast_cds[heg$gene_id])
#' rscu_heg <- est_rscu(cf_heg)
#' head(rscu_heg)
#'
est_rscu <- function(cf, weight = 1, pseudo_cnt = 1, codon_table = get_codon_table(),
                     level = 'subfam', incl_stop = FALSE){
    aa_code <- cts <- codon <- . <- rscu <- prop <- NULL # due to NSE notes in R CMD check
    if(!level %in% c('amino_acid', 'subfam')){
      stop('Possible values for `level` are "amino_acid" and "subfam"')
    }
    codon_freq <- colSums(cf * weight)
    if(!incl_stop){
      codon_table <- codon_table[aa_code != '*']
    }
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
#' @param domain The taxonomic domain of interest. "Eukarya" (default), "Bacteria" or "Archaea".
#' @param plot FALSE (default) or TRUE. Whether to keep the columns required for plotting.
#' @returns a data.table of codon-anticodon pairing information. The columns represent 
#' the pairing type, codon, corresponding anticodon, and the encoded amino acid when the argument "plot" is FALSE.
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
    # no pairing to stop codons or anticodon corresponding to stop codons
    ca_pairs <- ca_pairs[!codon %in% stop_codons]
    ca_pairs <- ca_pairs[!anticodon %in% as.character(rev_comp(stop_codons))]
    # only wobble among synonymous codons and anticodons
    ca_pairs[codon_table, codon_aa := i.aa_code, on = .(codon)]
    ca_pairs[codon_table, amino_acid := i.amino_acid, on = .(codon)]
    ca_pairs[codon_table, anticodon_aa := i.aa_code, on = .(anticodon)]
    ca_pairs <- ca_pairs[codon_aa == anticodon_aa]
    
    if(domain == "Bacteria"){
      dt <-data.table::data.table(codon_b1 = "A", codon_b2 = "T", type = "LA", base_codon = "A", 
                             base_anti = "G", codon = "ATA", anticodon = "CAT", codon_aa = "I", 
                             amino_acid = "Ile", anticodon_aa = "I")
      ca_pairs <- rbind(ca_pairs, dt)
    }else if(domain == "Archaea"){
      dt <-data.table::data.table(codon_b1 = "A", codon_b2 = "T", type = "agmA", base_codon = "A", 
                                  base_anti = "G", codon = "ATA", anticodon = "CAT", codon_aa = "I", 
                                  amino_acid = "Ile", anticodon_aa = "I")
      ca_pairs <- rbind(ca_pairs, dt)
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
#' @param pairs a table of codon-anticodon pairing derived from \code{ca_pairs}
#' @returns a plot on possible codon-anticodons pairings
#' @importFrom data.table ':='
#' @importFrom rlang .data
#' @export
#' @examples
#' # plot possible codon and anticodon pairings for the vertebrate mitochondrial genetic code
#' ctab <- get_codon_table(gcid = '2')
#' pairs <- ca_pairs(ctab, plot = TRUE)
#' plot_ca_pairs(ctab, pairs)
#'
#' # plot possible codon and anticodon pairings for the standard genetic code in bacteria
#' plot_ca_pairs(pairs = ca_pairs(domain = "Bacteria", plot = TRUE))
#' 
plot_ca_pairs <- function(codon_table = get_codon_table(), pairs = pairs){
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
      data = pairs,
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

#' Extract tRNA gene copy numbers from nature tRNA sequences
#'
#' \code{extract_trna_gcn} processes tRNA sequence data from  GtRNADB to 
#' extract gene copy numbers for each tRNA  type. This information is essential
#' for calculating tRNA availability weights used in TAI analysis.
#'
#' @param trna_seq A named vector or DNAStringSet of tRNA sequences, typically 
#'   from GtRNADB. Sequence names should follow the standard format containing 
#'   amino acid and anticodon information (e.g., "tRNA-Ala-AGC-1-1").
#' @return A named table of tRNA gene copy numbers. Names are in the format 
#'   "AminoAcid-Anticodon" (e.g., "Ala-AGC") and values represent the count 
#'   of genes encoding each tRNA type. Initiator tRNAs (iMet, fMet) and 
#'   undetermined tRNAs (Und-NNN) are automatically excluded as they serve 
#'   specialized functions in translation initiation.
#' @export
#' @examples
#' # Extract tRNA gene copy numbers for yeast
#' trna_gcn <- extract_trna_gcn(yeast_trna)
#' head(trna_gcn)
#' 
#' # View the distribution of tRNA gene copies
#' hist(trna_gcn, main = "Distribution of tRNA gene copy numbers")
#' 
extract_trna_gcn <- function(trna_seq){
  trna_copy <- sub('.*tRNA-(.*?)-\\d.*', '\\1', names(trna_seq))
  trna_copy <- trna_copy[!trna_copy %in% c("fMet-CAT", "iMet-CAT", "Und-NNN")]
  trna_gcn <- table(as.vector(trna_copy))
  return(trna_gcn)
}


#' Estimate tRNA weights for TAI calculation
#'
#' \code{est_trna_weight} calculates tRNA weights for each codon based on tRNA 
#' availability and codon-anticodon pairing efficiency. These weights are used 
#' in tRNA Adaptation Index (TAI) calculations and reflect how well each codon 
#' is supported by the cellular tRNA pool.
#'
#' @param trna_level A named numeric vector of tRNA expression levels or gene 
#'   copy numbers. Names should be in the format "AminoAcid-Anticodon" 
#'   (e.g., "Ala-GCA"). Each value represents the abundance of that tRNA species.
#' @param codon_table A codon table defining the genetic code, derived from 
#'   \code{get_codon_table()} or \code{create_codon_table()}.
#' @param domain Character string specifying the taxonomic domain: "Eukarya" 
#'   (default), "Bacteria", or "Archaea". This determines the codon-anticodon 
#'   pairing rules and selection penalties. Specify either "domain" or "s".
#' @param s A named list of selection penalties for non-Watson-Crick pairings. 
#'   If provided, overrides the default domain-specific penalties. Specify 
#'   either "domain" or "s".
#' @return A data.table containing comprehensive tRNA weight information with columns:
#'   \itemize{
#'     \item \code{aa_code}: Single-letter amino acid code
#'     \item \code{amino_acid}: Three-letter amino acid abbreviation
#'     \item \code{codon}: Codon sequence
#'     \item \code{subfam}: Codon subfamily identifier
#'     \item \code{anticodon}: Corresponding anticodon sequence
#'     \item \code{trna_id}: tRNA identifier (amino_acid-anticodon)
#'     \item \code{ac_level}: tRNA abundance level
#'     \item \code{W}: Absolute adaptiveness value
#'     \item \code{w}: Relative adaptiveness (normalized weight for TAI)
#'   }
#' @importFrom data.table ':='
#' @references dos Reis M, Savva R, Wernisch L. 2004. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res 32:5036-5044.
#' @references Sabi R, Tuller T. 2014. Modelling the efficiency of codon-tRNA interactions based on codon usage bias. DNA Res 21:511-526.
#' @export
#' @examples
#' # Calculate tRNA weights for yeast using gene copy numbers
#' yeast_trna_w <- est_trna_weight(yeast_trna_gcn)
#' head(yeast_trna_w)
#' 
#' # View the weight distribution
#' hist(yeast_trna_w$w, main = "Distribution of tRNA weights")
#' 
est_trna_weight <- function(trna_level, codon_table = get_codon_table(), domain = "Eukarya", s = NULL){
    anticodon <- aa_code <- ac_level <- penalty <- amino_acid <- NULL # due to NSE notes in R CMD check
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
    ca_pairs[s, penalty := i.values, on = .(type = ind)]
    ca_pairs[, trna_id := paste(amino_acid, anticodon, sep = "-")]
    ca_pairs <- ca_pairs[trna_id %in% names(trna_level)]
    ca_pairs[, ac_level := trna_level[trna_id]]

    dtt_W <- ca_pairs[, .(W = sum(ac_level * (1 - penalty))), by = .(codon)]
    codon_table[dtt_W, W := i.W, on = .(codon)]
    codon_table[, w := W/max(W, na.rm = TRUE)]

    # using geometric mean of w for rare cases that a codon has no matching tRNA
    # probably due to incomplete tRNA annotation
    w0 <- codon_table[!is.na(w), w]
    mean_w <- exp(mean(log(w0)))
    codon_table[is.na(w), w := mean_w]
    return(codon_table)
}


#' Identify optimal codons using statistical modeling
#'
#' \code{est_optimal_codons} identifies optimal codons within each codon family 
#' or amino acid group using binomial regression. Optimal codons are those whose 
#' usage correlates positively with high gene expression or negatively with 
#' codon usage bias (ENC), suggesting they are preferred for efficient translation.
#'
#' @param cf A matrix of codon frequencies as calculated by \code{count_codons()}.
#'   Rows represent sequences and columns represent codons.
#' @param codon_table A codon table defining the genetic code, derived from 
#'   \code{get_codon_table()} or \code{create_codon_table()}.
#' @param level Character string specifying the analysis level: "subfam" (default, 
#'   analyzes codon subfamilies) or "amino_acid" (analyzes at amino acid level).
#' @param gene_score A numeric vector of gene-level scores used to identify 
#'   optimal codons. Length must equal the number of rows in \code{cf}. Common 
#'   choices include:
#'   \itemize{
#'     \item Gene expression levels (RPKM, TPM, FPKM) - optionally log-transformed
#'     \item Protein abundance measurements
#'     \item Custom gene importance scores
#'   }
#'   If not provided, the negative of ENC values will be used (lower ENC = higher bias).
#' @param fdr Numeric value specifying the false discovery rate threshold for 
#'   determining statistical significance of codon optimality (default depends on method).
#' @return A data.table containing the input codon table with additional columns 
#'   indicating codon optimality status, statistical significance, and effect sizes 
#'   from the regression analysis. The columns include single-letter abbreviation 
#'   of the amino acid, three-letter abbreviation, codon, codon subfamily, 
#'   regression coefficient, regression P-value, Benjamini and Hochberg corrected 
#'   Q-value, and indication of whether the codon is optimal.
#' @references 
#' Presnyak V, Alhusaini N, Chen YH, Martin S, Morris N, Kline N, Olson S, 
#' Weinberg D, Baker KE, Graveley BR, et al. 2015. Codon optimality is a major 
#' determinant of mRNA stability. Cell 160:1111-1124.
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
#' @returns a data.table of codons and their CSCs. The columns include codon, codon stability coefficient, and correlation P-value.
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

#' Estimate Amino Acid Usage Frequencies of CDSs.
#'
#' @param cf matrix of codon frequencies as calculated by `count_codons()`.
#' @param codon_table codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @returns a data.table with amino acid frequencies of CDSs. The columns include three-letter abbreviation of the amino acid, 
#' single-letter abbreviation, usage frequency of the amino acid in all sequences, and usage frequency proportion.
#' @export
#' @examples
#' # estimate amino acid frequencies of yeast genes
#' cf_all <- count_codons(yeast_cds)
#' aau <- est_aau(cf_all)
#' print(aau)
est_aau <- function(cf, codon_table = get_codon_table()){
  aa_code <- count_codon <- count <- codon <- amino_acid <- aa_code <- proportion <- . <- NULL # due to NSE notes in R CMD check
  codon_table <- data.table::as.data.table(codon_table)
  codon_table <- codon_table[aa_code != '*']
  codon_table[, count_codon := colSums(cf)[codon]]
  aau <- codon_table[, .(count = sum(count_codon)), by = .(amino_acid, aa_code)]
  aau[, proportion := count/sum(aau$count)]
  return(aau)
}
