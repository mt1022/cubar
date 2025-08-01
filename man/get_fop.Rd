% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_index.R
\name{get_fop}
\alias{get_fop}
\title{Calculate fraction of optimal codons (Fop)}
\usage{
get_fop(cf, op = NULL, codon_table = get_codon_table(), ...)
}
\arguments{
\item{cf}{A matrix of codon frequencies as calculated by \code{count_codons()}.
Rows represent sequences and columns represent codons.}

\item{op}{A character vector specifying which codons are considered optimal.
If not provided, optimal codons will be determined automatically using
\code{est_optimal_codons()}.}

\item{codon_table}{A codon table defining the genetic code, derived from
\code{get_codon_table()} or \code{create_codon_table()}.}

\item{...}{Additional arguments passed to \code{est_optimal_codons()} when
optimal codons are determined automatically.}
}
\value{
A named numeric vector of Fop values (ranging from 0 to 1). Names
  correspond to sequence identifiers from the input matrix. Higher values
  indicate greater usage of optimal codons.
}
\description{
\code{get_fop} calculates the fraction of optimal codons (Fop) for each
coding sequence, which represents the proportion of codons that are
considered optimal for translation efficiency. Higher Fop values suggest
stronger selection for optimal codon usage.
}
\examples{
# Calculate Fop for yeast genes (optimal codons determined automatically)
cf_all <- count_codons(yeast_cds)
fop <- get_fop(cf_all)
head(fop)
hist(fop, main = "Distribution of Fop values")


}
\references{
Ikemura T. 1981. Correlation between the abundance of Escherichia coli transfer RNAs
  and the occurrence of the respective codons in its protein genes: a proposal for a synonymous
  codon choice that is optimal for the E. coli translational system. J Mol Biol 151:389-409.
}
