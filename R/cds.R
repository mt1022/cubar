seq_to_codons <- function(seq){
    codons <- IRanges::Views(seq, start = seq(1, length(seq), 3), width = 3)
    return(as.character(codons))
}

#' TODO: min CDS length and subseq operations
trim_cds <- function(cds, gcid = '1', min_len = 6, check_len = TRUE, check_start = TRUE,
                     check_stop = TRUE, check_istop = TRUE,
                     rm_start = TRUE, rm_stop = TRUE, rm_istop = TRUE){
    gct <- codon_table(gcid)
    stop_codons <- gct[aa_code == '*', codon]
    # if input is RNA sequences, convert to DNA
    if(class(cds) == 'RNAStringSet') cds <- Biostrings::DNAStringSet(cds)
    cds <- cds[width(cds) >= min_len]
    if(check_len){
        # CDS length is multiple of 3
        cds <- cds[(width(cds) %% 3) == 0]
    }
    if(check_start){
        cds <- cds[subseq(cds, 1, 3) == 'ATG']
    }
    if(check_stop){
        cds <- cds[subseq(cds, width(cds) - 2, width(cds)) %in% stop_codons]
    }
    if(check_istop){
        w_istop <- sapply(subseq(cds, 1, width(cds) - 3), function(x){
            any(seq_to_codons(x) %in% stop_codons)
        })
        cds <- cds[!w_istop]
    }
}
