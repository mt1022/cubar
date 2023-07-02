#' Convert CDS to codons
#'
#' \code{seq_to_codons} converts a coding sequence to a vector of codons
#'
#' @param seq DNAString, or an object that can be coerced to a DNAString
#' @returns a character vector of codons
seq_to_codons <- function(seq){
    if(class(seq) != 'DNAString'){
        seq <- DNAString(seq)
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
rev_comp <- function(seqs){
    if(class(seqs) != "DNAStringSet"){
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
#' @returns DNAStringSet of filtered (and trimmed) CDS sequences
check_cds <- function(seqs, codon_table = get_codon_table(), min_len = 6,
                      check_len = TRUE, check_start = TRUE, check_stop = TRUE,
                      check_istop = TRUE, rm_start = TRUE, rm_stop = TRUE,
                      start_codons = c("ATG")){
    stop_codons <- codon_table[aa_code == '*', codon]
    # if input is RNA sequences, convert to DNA
    if(class(seqs) == 'RNAStringSet'){
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
#' @param seqs CDS sequences, DNAStringSet
#' @returns matrix of codon (column) frequencies of each CDS (row)
count_codons <- function(seqs, ...){
    cf <- Biostrings::trinucleotideFrequency(seqs, step = 3, ...)
    rownames(cf) <- names(seqs)
    return(cf)
}
