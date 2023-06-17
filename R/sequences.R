#' Convert CDS to codons
#'
#' \code{seq_to_codons} convert a coding sequence to a vector of codons
#'
#' @param seq DNAString, or an object that can be coerced to a DNAString
#' @return a character vector of codons
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

#' Quality control of CDS
#'
#' \code{qc_cds} performs quality control of CDS sequences by filtering some
#' peculiar sequences and optionally remove start or stop codons.
#'
#' @param seqs input CDS sequences
#' @return filtered (and trimmed) CDS, DNAStringSet
qc_cds <- function(seqs, min_len = 6, check_len = TRUE, check_start = TRUE,
                     check_stop = TRUE, check_istop = TRUE,
                     rm_start = TRUE, rm_stop = TRUE, gcid = '1'){
    gct <- get_codon_table(gcid)
    stop_codons <- gct[aa_code == '*', codon]
    # if input is RNA sequences, convert to DNA
    if(class(seqs) == 'RNAStringSet') seqs <- Biostrings::DNAStringSet(seqs)
    seqs <- seqs[IRanges::width(seqs) >= min_len]
    # CDS length is multiple of 3?
    if(check_len){
        seqs <- seqs[(IRanges::width(seqs) %% 3) == 0]
    }
    # begin with start codon ATG?
    if(check_start){
        seqs <- seqs[Biostrings::subseq(seqs, 1, 3) == 'ATG']
    }
    # end with a stop codon
    if(check_stop){
        x <- IRanges::width(seqs)
        seqs <- seqs[Biostrings::subseq(seqs, x - 2, x) %in% stop_codons]
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

#' Count codon frequencies
#'
#' \code{count_codons} tabulate all the 64 codons in input CDSs
#'
#' @param seqs CDSs, DNAStringSet
#' @return matrix of codon (column) frequencies in each CDS (row)
count_codons <- function(seqs, ...){
    cf <- Biostrings::trinucleotideFrequency(seqs, step = 3, ...)
    rownames(cf) <- names(seqs)
    return(cf)
}
