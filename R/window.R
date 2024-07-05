# functions for sliding-window analysis

#' slide window interval generator
#'
#' \code{slide} generates a data.table with start, center, and end columns
#' for a sliding window analysis.
#'
#' @param from integer, the start of the sequence
#' @param to integer, the end of the sequence
#' @param step integer, the step size
#' @param before integer, the number of values before the center of a window
#' @param after integer, the number of values after the center of a window
#' @return data.table with start, center, and end columns
#' @examples
#' slide(1, 10, step = 2, before = 1, after = 1)
#'
slide <- function(from, to, step = 1, before = 0, after = 0){
    center <- seq(from + before, to - after, by = step)
    start <- center - before
    end <- center + after
    data.table(start = start, center = center, end = end)
}


#' sliding window of codons
#'
#' \code{slide_codon} generates a data.table with start, center, and end columns
#' for a sliding window analysis of codons.
#'
#' @param seq DNAString, the sequence
#' @param step integer, the step size
#' @param before integer, the number of codons before the center of a window
#' @param after integer, the number of codons after the center of a window
#' @return data.table with start, center, and end columns
#' @examples
#' x <- Biostrings::DNAString('ATCTACATAGCTACGTAGCTCGATGCTAGCATGCATCGTACGATCGTCGATCGTAG')
#' slide_codon(x, step = 3, before = 1, after = 1)
#'
slide_codon <- function(seq, step = 1, before = 0, after = 0){
    slen <- (length(seq) %/% 3) * 3
    slide(from = 1, to = slen, step = step*3,
          before = before*3, after = 2 + after*3)
}


#' apply a cub index to a sliding window
#'
#' \code{slide_apply} applies a function to a sliding window of codons.
#'
#' @param seq DNAString, the sequence
#' @param FUN function, the function to apply
#' @param step integer, the step size in number of codons
#' @param before integer, the number of codons before the center of a window
#' @param after integer, the number of codons after the center of a window
#' @return data.table with start, center, end, and codon usage index columns
#' @export
#' @examples
#' sw <- slide_apply(yeast_cds[[1]], get_enc, step = 1, before = 10, after = 10)
#' ggplot(sw, aes(x = ceiling(center/3), y = enc)) + geom_line() + geom_point()
#'
slide_apply <- function(seq, FUN, step = 1, before = 0, after = 0, ...){
    xw <- slide_codon(seq, step = step, before = before, after = after)
    xs <- Views(seq, start = xw$start, end = xw$end)
    xs <- as(xs, 'DNAStringSet')
    xw[, index := FUN(count_codons(xs), ...)]
    xw[]
}
