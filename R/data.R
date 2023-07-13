#' yeast CDS sequences
#'
#' CDSs of all protein-coding genes in Saccharomyces_cerevisiae
#'
#' @format a DNAStringSet of 6600 sequences
#' @source <https://ftp.ensembl.org/pub/release-107/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz>
#' @examples
#' head(yeast_cds)
"yeast_cds"


#' Half life of yeast mRNAs
#'
#' Half life of yeast mRNAs in Saccharomyces_cerevisiae calculated from rRNA-deleted total RNAs by Presnyak et al.
#'
#' @format a data.frame with 3888 rows and three columns:
#' \describe{
#' \item{gene_id}{gene id}
#' \item{gene_name}{gene name}
#' \item{half_life}{mRNA half life in minutes}
#' }
#' @source <https://doi.org/10.1016/j.cell.2015.02.029>
#' @references Presnyak V, Alhusaini N, Chen YH, et al. Codon optimality is a major determinant of mRNA stability. Cell. 2015;160(6):1111-1124. doi:10.1016/j.cell.2015.02.029
#' @examples
#' head(yeast_half_life)
"yeast_half_life"


#' yeast mRNA expression levels
#'
#' Yeast mRNA FPKM determined from rRNA-depleted (RiboZero) total RNA-Seq libraries. RUN1_0_WT and RUN2_0_WT (0 min after RNA Pol II repression) were averaged and used here.
#'
#' @format a data.frame with 6717 rows and three columns:
#' \describe{
#' \item{gene_id}{gene ID}
#' \item{gene_name}{gene name}
#' \item{fpkm}{mRNA expression level in Fragments per kilobase per million reads}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57385&format=file&file=GSE57385%5FRiboZero%5FFPKM%2Etxt%2Egz>
#' @references Presnyak V, Alhusaini N, Chen YH, et al. Codon optimality is a major determinant of mRNA stability. Cell. 2015;160(6):1111-1124. doi:10.1016/j.cell.2015.02.029
#' @examples
#' head(yeast_exp)
"yeast_exp"


#' yeast tRNA gene copy numbers (GCN)
#'
#' Yeast tRNA gene copy numbers (GCN) by anticodon obtained from gtRNAdb.
#'
#' @format a named vector with a length of 41. Value names are anticodons.
#' @source <http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/sacCer3-mature-tRNAs.fa>
#' @references Chan, P.P. & Lowe, T.M. (2016) GtRNAdb 2.0: an expanded database of transfer RNA genes identified in complete and draft genomes. Nucl. Acids Res. 44(Database issue):D184-D189.Chan, P.P. & Lowe, T.M. (2016) GtRNAdb 2.0: an expanded database of transfer RNA genes identified in complete and draft genomes. Nucl. Acids Res. 44(Database issue):D184-D189.
#' @examples
#' yeast_trna_gcn
"yeast_trna_gcn"


#' amino acids to codons
#'
#'  A data.frame of mapping from amino acids to codons
#'
#' @format a data.frame with two columns: amino_acid, and codon.
#' \describe{
#' \item{amino_acid}{amino acid corresponding to the codon}
#' \item{codon}{codon identity}
#' }
#' @source It is actually the standard genetic code.
#' @examples
#' aa2codon
#'
"aa2codon"
