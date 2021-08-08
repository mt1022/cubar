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
#' \code{qc_cds} perform quality control of CDS sequences by filtering some peculiar sequences
#' and optionally remove start or stop codons.
#'
#' @param cds input CDS sequences
#' @return filtered (and trimmed) CDS, DNAStringSet
qc_cds <- function(cds, min_len = 6, check_len = TRUE, check_start = TRUE,
                     check_stop = TRUE, check_istop = TRUE,
                     rm_start = TRUE, rm_stop = TRUE, gcid = '1'){
    gct <- codon_table(gcid)
    stop_codons <- gct[aa_code == '*', codon]
    # if input is RNA sequences, convert to DNA
    if(class(cds) == 'RNAStringSet') cds <- Biostrings::DNAStringSet(cds)
    cds <- cds[IRanges::width(cds) >= min_len]
    # CDS length is multiple of 3?
    if(check_len){
        cds <- cds[(IRanges::width(cds) %% 3) == 0]
    }
    # begin with start codon ATG?
    if(check_start){
        cds <- cds[Biostrings::subseq(cds, 1, 3) == 'ATG']
    }
    # end with a stop codon
    if(check_stop){
        cds <- cds[Biostrings::subseq(cds, IRanges::width(cds) - 2, IRanges::width(cds)) %in% stop_codons]
    }
    # no internal stop codons?
    if(check_istop){
        w_istop <- sapply(Biostrings::subseq(cds, 1, IRanges::width(cds) - 3), function(x){
            any(seq_to_codons(x) %in% stop_codons)
        })
        cds <- cds[!w_istop]
    }
    # trimming
    if(rm_start){
        cds <- Biostrings::subseq(cds, 4, IRanges::width(cds))
    }
    if(rm_stop){
        cds <- Biostrings::subseq(cds, 1, IRanges::width(cds) - 3)
    }
    return(cds)
}
