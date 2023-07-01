#' Reverse complement
#'
#' \code{rev_comp} creates reverse complemented version of the input sequence
#'
#' @param seqs input sequences, DNAStringSet or named vector of sequences
#' @returns reverse complemented input sequences as a DNAStringSet.
rev_comp <- function(seqs){
    if(class(seqs) != "DNAStringSet"){
        seqs <- Biostrings::DNAStringSet(seqs)
    }
    Biostrings::reverseComplement(seqs)
}


#' Estimate RSCU
#'
#' \code{est_rscu} returns the RSCU value of codons
#'
#' @param seqs CDS sequences
#' @param weight a vector of the same length as `seqs` that gives different weights to
#'   CDSs when count codons. for example, it could be gene expression levels.
#' @param pseudo_cnt pseudo count to avoid dividing by zero. This may occur when
#'   only a few sequences are available for RSCU calculation.
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`
#' @returns a data.table of codon info and RSCU values
#' @references
est_rscu <- function(seqs, weight = 1, pseudo_cnt = 1, codon_table = get_codon_table()){
    seqs <- Biostrings::DNAStringSet(seqs)
    codon_freq <- colSums(count_codons(seqs) * weight)

    codon_table <- codon_table[aa_code != '*']
    codon_table[, cts := codon_freq[codon]]
    codon_table[, `:=`(
        rscu = (cts + pseudo_cnt) / sum(cts + pseudo_cnt),
        w_cai = (cts + pseudo_cnt) / max(cts + pseudo_cnt)), by = .(subfam)]
    return(codon_table[])
}


#' Plot codon-anticodon pairing relationship
#'
#' \code{plot_ca_pairing} returns the RSCU value of codons
#'
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @param plot whether to plot the pairing relationship
#' @importFrom rlang .data
#' @returns a data.table of codon info and RSCU values
plot_ca_pairing <- function(codon_table = get_codon_table(), plot = TRUE){
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
    # no wobble for start codon
    ca_pairs <- ca_pairs[!(codon == 'ATG' & anticodon == 'TAT')]

    if(plot){
        ggplot2::ggplot(codon_table, ggplot2::aes(y = .data$codon_b3)) +
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
#' TODO test
#'
#' @param trna_level, named vector of tRNA level (or gene copy numbers), one value for each anticodon.
#'   vector names are anticodons.
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @param s list of non-Waston-Crick pairing panelty.
#' @return data.table of tRNA expression information
est_trna_weight <- function(trna_level, codon_table = get_codon_table(),
                            s = list(WC=0, IU=0, IC=0.4659, IA=0.9075, GU=0.7861, UG=0.6295)){
    codon_table[, anticodon := as.character(Biostrings::reverseComplement(
        Biostrings::DNAStringSet(codon_table$codon)))]
    codon_table <- codon_table[aa_code != '*']

    codon_table[, ac_level := trna_level[anticodon]]
    codon_table[is.na(ac_level), ac_level := 0]

    ca_pairs <- show_ca_pairs(plot = FALSE)
    s <- stack(s)
    ca_pairs[s, penality := i.values, on = .(type = ind)]
    ca_pairs <- ca_pairs[anticodon %in% names(trna_level)]
    ca_pairs[, ac_level := trna_level[anticodon]]
    dtt_W <- ca_pairs[, .(W = sum(ac_level * (1 - penality))), by = .(codon)]
    codon_table[dtt_W, W := i.W, on = .(codon)]
    codon_table[, w := W/max(W)]
    mean_w <- mean(codon_table$w, na.rm = TRUE)
    codon_table[is.na(w), w := mean_w]
    return(codon_table)
}


#' Estimate optimal codons
#'
#' \code{est_toptimal_codons} determine optimal codon of each codon family with binomial regression.
#'
#' @param seqs CDS sequences of all protein-coding genes. One for each gene.
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @returns data.table of optimal codons
est_optimal_codons <- function(seqs, ctab = get_codon_table()){
    enc <- get_enc(seqs, codon_table = ctab)
    cf_all <- count_codons(seqs)
    # regression analysis for each codon sub-family
    binreg <- lapply(split(ctab$codon, f = ctab$subfam), function(x){
        cf <- cf_all[, x, drop = FALSE]
        if(ncol(cf) == 1){
            data.table::data.table(
                codon = colnames(cf), coef = 0, se = 0, zvalue = 0, pvalue = 0)
        }else{
            total <- rowSums(cf)
            res <- apply(cf, 2, function(x){
                x
                fit <- glm(cbind(x, total - x) ~ enc, family = 'binomial')
                summary(fit)$coefficients[-1, ]
            })
            res <- data.table::as.data.table(
                as.data.frame(t(res)), keep.rownames = 'codon')
            data.table::setnames(res, c('codon', 'coef', 'se', 'zvalue', 'pvalue'))
        }
    })
    bingreg <- data.table::rbindlist(binreg, idcol = 'subfam')
    bingreg <- ctab[bingreg, on = .(codon, subfam)]
    return(bingreg)
}

#' Estimate Codon Stabilization Coefficient
#'
#' \code{get_csc} calculate codon occurrence to mRNA stability correlation coefficients (Default to Pearson's).
#'
#' @param seqs CDS sequences of all protein-coding genes. One for each gene.
#' @param half_life data.frame of mRNA half life (gene_id & half_life are column names).
#' @param codon_table a table of genetic code derived from `get_codon_table` or `create_codon_table`.
#' @param cor_method method name passed to `cor.test` used for calculating correlation coefficients.
#' @returns data.table of optimal codons.
est_csc <- function(seqs, half_life, codon_table = get_codon_table(), cor_method = 'pearson'){
    non_stop_codons <- codon_table[aa_code != '*', codon]
    cf <- count_codons(seqs)

    common_genes <- intersect(half_life$gene_id, rownames(cf))
    half_life <- half_life[half_life$gene_id %in% common_genes, ]
    cf <- cf[half_life$gene_id, non_stop_codons]
    cp <- cf / rowSums(cf)  # codon proportions

    csc <- apply(cp, 2, function(prop){
        res <- cor.test(prop, half_life$half_life, method = cor_method)
        c(corr = res$estimate, pvalue = res$p.value)
    })
    csc <- data.table::data.table(t(csc), keep.rownames = TRUE)
    data.table::setnames(csc, c('codon', 'csc', 'pvalue'))
    return(csc)
}
