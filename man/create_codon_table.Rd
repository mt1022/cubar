% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codon_info.R
\name{create_codon_table}
\alias{create_codon_table}
\title{create custom codon table from a data frame}
\usage{
create_codon_table(aa2codon)
}
\arguments{
\item{aa2codon}{a data frame with two columns: amino_acid (Ala, Arg, etc.) and codon.}
}
\value{
a \code{data.table} with four columns: aa_code, amino_acid, codon, and subfam.
}
\description{
\code{create_codon_table} creates codon table from data frame of aa to codon mapping.
}
\examples{
head(aa2codon)
create_codon_table(aa2codon = aa2codon)
}
