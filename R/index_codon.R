#' Reverse complement
rev_comp <- function(seqs){
    Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs))
}

#' Calculate RSCU
#'
#' \code{get_rscu} returns the RSCU value of codons
#'
#' @param seqs CDS sequences
#' @weight a vector of the same length as \code{seqs} that gives different
#' weights to CDSs when count codons. It could be gene expression levels.
#' @pseudo_cnt pseudo count to avoid dividing by zero. This may occur when
#' only a few sequences are available for RSCU calculation.
#' @gcid ID or name of genetic code. Support for non-standard genetic code will
#' be added in the future.
#' @return a data.table
#' @importFrom rlang .data
#' @references
get_rscu <- function(seqs, weight = 1, pseudo_cnt = 1, gcid = '1'){
    seqs <- Biostrings::DNAStringSet(seqs)
    codon_freq <- colSums(count_codons(seqs) * weight)

    codon_table <- get_codon_table(gcid)
    codon_table <- codon_table[aa_code != '*']
    codon_table[, cts := codon_freq[codon]]
    codon_table[, `:=`(
        rscu = (cts + pseudo_cnt) / sum(cts + pseudo_cnt),
        w_cai = (cts + pseudo_cnt) / max(cts + pseudo_cnt)), by = .(subfam)]
    return(codon_table[])
}

show_ca_pairs <- function(gcid='1', plot = TRUE){
    codon_table <- get_codon_table(gcid)
    codon_table[, anticodon := as.character(rev_comp(codon_table$codon))]
    codon_table[, c('codon_b1', 'codon_b2', 'codon_b3') :=
                    data.table::tstrsplit(codon, '')]
    bases <- c('T', 'C', 'A', 'G')
    codon_table[, codon_b1 := factor(codon_b1, levels = bases)]
    codon_table[, ]
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
        # note: base_anti here is used for visualization, not referring to 1st
        #   base of anticodon
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
    ca_pairs <- ca_pairs[!codon %in% stop_codons]
    ca_pairs <- ca_pairs[!anticodon %in% as.character(rev_comp(stop_codons))]
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

#' Calculate tRNA w
#'
#' \code{get_trna_weight} compute the tRNA weight per codon for TAI calculation.
#' This weight reflects relative tRNA availability for each codon.
#'
#' TODO test
#'
#' @param trna_level, named vector of tRNA level, one value for each anticodon.
#'   vector names are anticodons.
#' @return data.table of tRNA expression information
get_trna_weight <- function(trna_level, gcid = '1', s = list(
    WC=0, IU=0, IC=0.4659, IA=0.9075, GU=0.7861, UG=0.6295)){
    codon_table <- get_codon_table(gcid)
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
