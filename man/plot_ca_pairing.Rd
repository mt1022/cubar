% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codon_weight.R
\name{plot_ca_pairing}
\alias{plot_ca_pairing}
\title{Plot codon-anticodon pairing relationship}
\usage{
plot_ca_pairing(codon_table = get_codon_table(), plot = TRUE)
}
\arguments{
\item{codon_table}{a table of genetic code derived from \code{get_codon_table} or
\code{create_codon_table}.}

\item{plot}{whether to plot the pairing relationship}
}
\value{
a data.table of codon info and RSCU values
}
\description{
\code{plot_ca_pairing} show possible codon-anticodons pairings
}
\examples{
ctab <- get_codon_table(gcid = '2')
pairing <- plot_ca_pairing(ctab)
head(pairing)

}
