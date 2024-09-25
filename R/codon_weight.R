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


#' Plot codon-anticodon pairing relationship
#'
#' \code{plot_ca_pairing} show possible codon-anticodons pairings
#'
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or
#'   \code{create_codon_table}.
#' @param plot whether to plot the pairing relationship
#' @returns a data.table of codon info and RSCU values
#' @importFrom data.table ':='
#' @importFrom rlang .data
#' @export
#' @examples
#' ctab <- get_codon_table(gcid = '2')
#' pairing <- plot_ca_pairing(ctab)
#' head(pairing)
#'
plot_ca_pairing <- function(codon_table = get_codon_table(), plot = TRUE){
    anticodon <- codon <- codon_b1 <- codon_b2 <- codon_b3 <- NULL # due to NSE notes in R CMD check
    . <- aa_code <- base_codon <- base_anti <- type <- NULL
    anticodon_aa <- codon_aa <- i.aa_code <- NULL
    codon_table <- data.table::copy(codon_table)
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
    ca_pairs[codon_table, anticodon_aa := i.aa_code, on = .(anticodon)]
    ca_pairs <- ca_pairs[codon_aa == anticodon_aa]

    if(plot){
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
        # return silently
        invisible(ca_pairs[, .(type, codon, anticodon)])
    }else{
        return(ca_pairs[, .(type, codon, anticodon)])
    }
}


#' Estimate tRNA weight w
#'
#' \code{est_trna_weight} compute the tRNA weight per codon for TAI calculation.
#' This weight reflects relative tRNA availability for each codon.
#'
#' @param trna_level, named vector of tRNA level (or gene copy numbers), one value for each anticodon.
#'   vector names are anticodons.
#' @param codon_table a table of genetic code derived from \code{get_codon_table} or \code{create_codon_table}.
#' @param s list of non-Waston-Crick pairing panelty.
#' @returns data.table of tRNA expression information.
#' @importFrom data.table ':='
#' @references dos Reis M, Savva R, Wernisch L. 2004. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res 32:5036-5044.
#' @export
#' @examples
#' # estimate codon tRNA weight for yeasts
#' est_trna_weight(yeast_trna_gcn)
#'
est_trna_weight <- function(trna_level, codon_table = get_codon_table(),
                            s = list(WC=0, IU=0, IC=0.4659, IA=0.9075, GU=0.7861, UG=0.6295)){
    anticodon <- aa_code <- ac_level <- penality <- NULL # due to NSE notes in R CMD check
    i.values <- . <- ind <- codon <- W <- i.W <- w <- NULL # due to NSE notes in R CMD check
    codon_table[, anticodon := as.character(Biostrings::reverseComplement(
        Biostrings::DNAStringSet(codon_table$codon)))]
    codon_table <- codon_table[aa_code != '*']

    codon_table[, ac_level := trna_level[anticodon]]
    codon_table[is.na(ac_level), ac_level := 0]

    ca_pairs <- plot_ca_pairing(codon_table = codon_table, plot = FALSE)
    s <- utils::stack(s)
    ca_pairs[s, penality := i.values, on = .(type = ind)]
    ca_pairs <- ca_pairs[anticodon %in% names(trna_level)]
    ca_pairs[, ac_level := trna_level[anticodon]]
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
#' \code{est_toptimal_codons} determine optimal codon of each codon family with binomial regression.
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
