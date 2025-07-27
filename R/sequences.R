#' Convert a coding sequence to a codon vector
#'
#' \code{seq_to_codons} converts a coding sequence (CDS) into a vector of codons by 
#' splitting the sequence into non-overlapping triplets starting from the first position.
#'
#' @param seq A coding sequence as a DNAString object, or any object that can be 
#'   coerced to a DNAString (e.g., character string).
#' @return A character vector where each element represents a codon (3-nucleotide sequence).
#' @export
#' @examples
#' # Convert a CDS sequence to a sequence of codons
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


#' Generate reverse complement sequences
#'
#' \code{rev_comp} generates the reverse complement of input DNA sequences. 
#' This is commonly used for analyzing complementary strands or anticodon sequences.
#'
#' @param seqs Input DNA sequences as a DNAStringSet object, or a named vector 
#'   of sequences that can be coerced to DNAStringSet.
#' @return A DNAStringSet object containing the reverse complemented sequences.
#' @export
#' @examples
#' # Reverse complement of codons
#' rev_comp(Biostrings::DNAStringSet(c('TAA', 'TAG')))
#'
rev_comp <- function(seqs){
    if(!inherits(seqs, "DNAStringSet")){
        seqs <- Biostrings::DNAStringSet(seqs)
    }
    Biostrings::reverseComplement(seqs)
}


#' Quality control and preprocessing of coding sequences
#'
#' \code{check_cds} performs comprehensive quality control on coding sequences (CDS) 
#' by filtering sequences based on various criteria and optionally removing start 
#' or stop codons. This function ensures that sequences meet the requirements for 
#' downstream codon usage analysis.
#'
#' @param seqs Input CDS sequences as a DNAStringSet or compatible object.
#' @param min_len Minimum CDS length in nucleotides (default: 6).
#' @param check_len Logical. Check whether CDS length is divisible by 3 (default: TRUE).
#' @param check_start Logical. Check whether CDSs begin with valid start codons (default: TRUE).
#' @param check_stop Logical. Check whether CDSs end with valid stop codons (default: TRUE).
#' @param check_istop Logical. Check for internal stop codons (default: TRUE).
#' @param rm_start Logical. Remove start codons from the sequences (default: TRUE).
#' @param rm_stop Logical. Remove stop codons from the sequences (default: TRUE).
#' @param codon_table Codon table matching the genetic code of the input sequences.
#'   Generated using \code{get_codon_table()} or \code{create_codon_table()}.
#' @param start_codons Character vector specifying valid start codons (default: "ATG").
#' @return A DNAStringSet containing filtered and optionally trimmed CDS sequences 
#'   that pass all quality control checks.
#' @export
#' @examples
#' # Perform CDS sequence quality control for a sample of yeast genes
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


#' Count codon frequencies in coding sequences
#'
#' \code{count_codons} tabulates the frequency of all 64 possible codons across 
#' input coding sequences. This function provides the foundation for most codon 
#' usage bias analyses in the cubar package.
#'
#' @param seqs Coding sequences as a DNAStringSet object, or compatible input 
#'   that can be coerced to DNAStringSet.
#' @param ... Additional arguments passed to \code{Biostrings::trinucleotideFrequency}.
#' @return A matrix where rows represent individual CDS sequences and columns 
#'   represent the 64 possible codons. Each cell contains the frequency count 
#'   of the corresponding codon in the respective sequence.
#' @export
#' @examples
#' # Count codon frequencies across all yeast CDS sequences
#' cf_all <- count_codons(yeast_cds)
#' dim(cf_all)
#' cf_all[1:5, 1:5]
#' 
#' # Count codons for a single sequence
#' count_codons(yeast_cds[1])
#'
count_codons <- function(seqs, ...){
    if(!inherits(seqs, 'DNAStringSet')){
        seqs <- Biostrings::DNAStringSet(seqs)
    }
    cf <- Biostrings::trinucleotideFrequency(seqs, step = 3, ...)
    rownames(cf) <- names(seqs)
    return(cf)
}
