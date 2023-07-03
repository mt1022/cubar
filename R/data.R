#' yeast CDS sequences
#'
#' CDSs of all protein-coding genes in Saccharomyces_cerevisiae
#'
#' @format ## `yeast_cds`
#' A DNAStringSet of 6600 yeast CDS sequences
#' @source <https://ftp.ensembl.org/pub/release-107/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz>
"yeast_cds"

#' Half life of yeast mRNAs
#'
#' Half life of yeast mRNAs in Saccharomyces_cerevisiae calculated from rRNA-deleted total RNAs by Presnyak et al.
#'
#' @format ## `yeast_half_life`
#' A data.frame with 3888 rows and three columns (gene_id, gene_name, half_life in min)
#' @source <https://doi.org/10.1016/j.cell.2015.02.029>
"yeast_half_life"


#' amino acids to codons
#'
#'  A data.frame of mapping from amino acids to codons
#'
#' @format ## `aa2codon`
#' data.frame with two columns: amino_acid, and codons.
#' @source It is actually the standard genetic code.
"aa2codon"
