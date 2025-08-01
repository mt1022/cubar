% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codon_info.R
\name{get_codon_table}
\alias{get_codon_table}
\title{Retrieve codon table by NCBI genetic code ID}
\usage{
get_codon_table(gcid = "1")
}
\arguments{
\item{gcid}{A character string specifying the NCBI genetic code ID. Use 
\code{show_codon_tables()} to view all available genetic codes and their 
corresponding IDs. Default is "1" (standard genetic code).}
}
\value{
A data.table with four columns:
  \itemize{
    \item \code{aa_code}: Single-letter amino acid code
    \item \code{amino_acid}: Three-letter amino acid abbreviation  
    \item \code{codon}: Three-nucleotide codon sequence
    \item \code{subfam}: Codon subfamily identifier (amino_acid_XY format)
  }
}
\description{
\code{get_codon_table} creates a standardized codon table based on genetic 
codes cataloged by NCBI. This function provides the mapping between codons 
and amino acids for different organisms and organelles, which is essential 
for accurate codon usage analysis.
}
\examples{
# Standard genetic code (used by most organisms)
standard_code <- get_codon_table()
head(standard_code)

# Vertebrate mitochondrial genetic code
mito_code <- get_codon_table(gcid = '2')
head(mito_code)

}
