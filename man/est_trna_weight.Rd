% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codon_weight.R
\name{est_trna_weight}
\alias{est_trna_weight}
\title{Estimate tRNA weights for TAI calculation}
\usage{
est_trna_weight(
  trna_level,
  codon_table = get_codon_table(),
  domain = "Eukarya",
  s = NULL
)
}
\arguments{
\item{trna_level}{A named numeric vector of tRNA expression levels or gene 
copy numbers. Names should be in the format "AminoAcid-Anticodon" 
(e.g., "Ala-GCA"). Each value represents the abundance of that tRNA species.}

\item{codon_table}{A codon table defining the genetic code, derived from 
\code{get_codon_table()} or \code{create_codon_table()}.}

\item{domain}{Character string specifying the taxonomic domain: "Eukarya" 
(default), "Bacteria", or "Archaea". This determines the codon-anticodon 
pairing rules and selection penalties. Specify either "domain" or "s".}

\item{s}{A named list of selection penalties for non-Watson-Crick pairings. 
If provided, overrides the default domain-specific penalties. Specify 
either "domain" or "s".}
}
\value{
A data.table containing comprehensive tRNA weight information with columns:
  \itemize{
    \item \code{aa_code}: Single-letter amino acid code
    \item \code{amino_acid}: Three-letter amino acid abbreviation
    \item \code{codon}: Codon sequence
    \item \code{subfam}: Codon subfamily identifier
    \item \code{anticodon}: Corresponding anticodon sequence
    \item \code{trna_id}: tRNA identifier (amino_acid-anticodon)
    \item \code{ac_level}: tRNA abundance level
    \item \code{W}: Absolute adaptiveness value
    \item \code{w}: Relative adaptiveness (normalized weight for TAI)
  }
}
\description{
\code{est_trna_weight} calculates tRNA weights for each codon based on tRNA 
availability and codon-anticodon pairing efficiency. These weights are used 
in tRNA Adaptation Index (TAI) calculations and reflect how well each codon 
is supported by the cellular tRNA pool.
}
\examples{
# Calculate tRNA weights for yeast using gene copy numbers
yeast_trna_w <- est_trna_weight(yeast_trna_gcn)
head(yeast_trna_w)

# View the weight distribution
hist(yeast_trna_w$w, main = "Distribution of tRNA weights")

}
\references{
dos Reis M, Savva R, Wernisch L. 2004. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res 32:5036-5044.

Sabi R, Tuller T. 2014. Modelling the efficiency of codon-tRNA interactions based on codon usage bias. DNA Res 21:511-526.
}
