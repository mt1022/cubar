% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_index.R
\name{get_gc}
\alias{get_gc}
\title{Calculate GC content of coding sequences}
\usage{
get_gc(cf)
}
\arguments{
\item{cf}{A matrix of codon frequencies as calculated by \code{count_codons()}.
Rows represent sequences and columns represent codons.}
}
\value{
A named numeric vector of GC content values (ranging from 0 to 1).
  Names correspond to sequence identifiers from the input matrix.
}
\description{
\code{get_gc} calculates the overall GC content (percentage of guanine and
cytosine nucleotides) for each coding sequence.
}
\examples{
# Calculate GC content for yeast genes
cf_all <- count_codons(yeast_cds)
gc <- get_gc(cf_all)
head(gc)
hist(gc, main = "Distribution of GC content")

}
