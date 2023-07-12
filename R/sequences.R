#' Convert CDS to codons
#'
#' \code{seq_to_codons} converts a coding sequence to a vector of codons
#'
#' @param seq DNAString, or an object that can be coerced to a DNAString
#' @returns a character vector of codons
#' @export
#' @examples
#' # convert a CDS sequence to a sequence of codons
#' seq_to_codons('ATGTGGTAG')
#' seq_to_codons(yeast_cds[[1]])
#'
seq_to_codons <- function(seq){
    if(!inherits(seq, 'DNAString')){
        seq <- Biostrings::DNAString(seq)
    }
    if(length(seq) < 3){
        stop('Input sequence too short')
    }
    codons <- IRanges::Views(seq, start = seq(1, length(seq) - 2, 3), width = 3)
    return(as.character(codons))
}


#' Reverse complement
#'
#' \code{rev_comp} creates reverse complemented version of the input sequence
#'
#' @param seqs input sequences, DNAStringSet or named vector of sequences
#' @returns reverse complemented input sequences as a DNAStringSet.
#' @export
#' @examples
#' # reverse complement of codons
#' rev_comp(Biostrings::DNAStringSet(c('TAA', 'TAG')))
#'
rev_comp <- function(seqs){
    if(!inherits(seqs, "DNAStringSet")){
        seqs <- Biostrings::DNAStringSet(seqs)
    }
    Biostrings::reverseComplement(seqs)
}


#' Quality control of CDS
#'
#' \code{check_cds} performs quality control of CDS sequences by filtering some
#' peculiar sequences and optionally remove start or stop codons.
#'
#' @param seqs input CDS sequences
#' @param min_len minimum CDS length in nt
#' @param check_len check whether CDS length is divisible by 3
#' @param check_start check whether CDSs have start codons
#' @param check_stop check whether CDSs have stop codons
#' @param check_istop check internal stop codons
#' @param rm_start whether to remove start codons
#' @param rm_stop whether to remove stop codons
#' @param codon_table codon table matching the genetic code of \code{seqs}
#' @param start_codons vector of start codons
#' @returns DNAStringSet of filtered (and trimmed) CDS sequences
#' @export
#' @examples
#' # CDS sequence QC for a sample of yeast genes
#' s <- head(yeast_cds, 10)
#' print(s)
#' check_cds(s)
#'
check_cds <- function(seqs, codon_table = get_codon_table(), min_len = 6,
                      check_len = TRUE, check_start = TRUE, check_stop = TRUE,
                      check_istop = TRUE, rm_start = TRUE, rm_stop = TRUE,
                      start_codons = c("ATG")){
    aa_code <- codon <- NULL  # due to NSE notes in R CMD check
    stop_codons <- codon_table[aa_code == '*', codon]
    # if input is RNA sequences, convert to DNA
    if(inherits(seqs, 'RNAStringSet')){
        seqs <- Biostrings::DNAStringSet(seqs)
    }
    seqs <- seqs[IRanges::width(seqs) >= min_len]
    # CDS length is divisible by 3?
    if(check_len){
        seqs <- seqs[(IRanges::width(seqs) %% 3) == 0]
    }
    # begin with start codon ATG?
    if(check_start){
        seqs <- seqs[as.character(Biostrings::subseq(seqs, 1, 3)) %in% start_codons]
    }
    # end with a stop codon
    if(check_stop){
        x <- IRanges::width(seqs)
        seqs <- seqs[as.character(Biostrings::subseq(seqs, x - 2, x)) %in% stop_codons]
    }
    # no internal stop codons?
    if(check_istop){
        x <- IRanges::width(seqs) - 3
        w_istop <- sapply(
            Biostrings::subseq(seqs, 1, x),
            function(seq){ any(seq_to_codons(seq) %in% stop_codons) })
        seqs <- seqs[!w_istop]
    }
    # trimming
    if(rm_start){
        seqs <- Biostrings::subseq(seqs, 4, IRanges::width(seqs))
    }
    if(rm_stop){
        seqs <- Biostrings::subseq(seqs, 1, IRanges::width(seqs) - 3)
    }
    return(seqs)
}


#' Count occurrences of different codons
#'
#' \code{count_codons} tabulates the occurrences of all the 64 codons in input CDSs
#'
#' @param seqs CDS sequences, DNAStringSet.
#' @param ... additional arguments passed to `Biostrings::trinucleotideFrequency`.
#' @returns matrix of codon (column) frequencies of each CDS (row).
#' @export
#' @examples
#' # count codon occurrences
#' cf_all <- count_codons(yeast_cds)
#' dim(cf_all)
#' cf_all[1:5, 1:5]
#' count_codons(yeast_cds[1])
#'
count_codons <- function(seqs, ...){
    cf <- Biostrings::trinucleotideFrequency(seqs, step = 3, ...)
    rownames(cf) <- names(seqs)
    return(cf)
}
