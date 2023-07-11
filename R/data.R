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
#' @references Presnyak V, Alhusaini N, Chen YH, et al. Codon optimality is a major determinant of mRNA stability. Cell. 2015;160(6):1111-1124. doi:10.1016/j.cell.2015.02.029
"yeast_half_life"


#' yeast mRNA expression levels
#'
#' Yeast mRNA FPKM determined from rRNA-depleted (RiboZero) total RNA-Seq libraries. RUN1_0_WT and RUN2_0_WT (0 min after RNA Pol II repression) were averaged and used here.
#'
#' @format ## `yeast_exp`
#' A data.frame with 6717 rows and three columns (gene_id, gene_name, fpkm)
#' @source <https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57385&format=file&file=GSE57385%5FRiboZero%5FFPKM%2Etxt%2Egz>
#' @references Presnyak V, Alhusaini N, Chen YH, et al. Codon optimality is a major determinant of mRNA stability. Cell. 2015;160(6):1111-1124. doi:10.1016/j.cell.2015.02.029
"yeast_exp"


#' amino acids to codons
#'
#'  A data.frame of mapping from amino acids to codons
#'
#' @format ## `aa2codon`
#' data.frame with two columns: amino_acid, and codons.
#' @source It is actually the standard genetic code.
"aa2codon"
