% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{codon_optimize}
\alias{codon_optimize}
\title{Optimize codon usage in coding sequences}
\usage{
codon_optimize(
  seq,
  optimal_codons = optimal_codons,
  cf = NULL,
  codon_table = get_codon_table(),
  level = "subfam",
  method = "naive",
  num_sequences = 1,
  organism = NULL,
  envname = "cubar_env",
  attention_type = "original_full",
  deterministic = TRUE,
  temperature = 0.2,
  top_p = 0.95,
  match_protein = FALSE,
  spliceai = FALSE
)
}
\arguments{
\item{seq}{A coding sequence as a DNAString object or any object that can 
be coerced to DNAString. The sequence should not include stop codons.}

\item{optimal_codons}{A table of optimal codons as generated by 
\code{est_optimal_codons()}, containing optimality information for each codon.}

\item{cf}{Matrix of codon frequencies from \code{count_codons()}. Required 
for "IDT" method to determine codon frequency distributions.}

\item{codon_table}{A codon table defining the genetic code, derived from 
\code{get_codon_table()} or \code{create_codon_table()}.}

\item{level}{Character string specifying optimization level: "subfam" (default, 
within codon subfamilies) or "amino_acid" (within amino acid groups). 
Required for "naive" and "IDT" methods.}

\item{method}{Character string specifying the optimization algorithm:
\itemize{
  \item \code{"naive"} (default): Simple replacement with optimal codons
  \item \code{"IDT"}: Method from Integrated DNA Technologies tool
  \item \code{"CodonTransformer"}: Neural network-based optimization
}}

\item{num_sequences}{Integer. Number of different optimized sequences to 
generate (default: 1). For "CodonTransformer" with deterministic=FALSE, 
each sequence is independently sampled.}

\item{organism}{Organism identifier (integer ID or string name) for 
"CodonTransformer" method. Must be from ORGANISM2ID in CodonUtils 
(e.g., "Escherichia coli general").}

\item{envname}{Environment name for "CodonTransformer" method. Should match 
your conda environment name (default: "cubar_env").}

\item{attention_type}{Attention mechanism type for "CodonTransformer":
\itemize{
  \item \code{"original_full"} (default): Standard attention
  \item \code{"block_sparse"}: Memory-efficient sparse attention
}}

\item{deterministic}{Logical. For "CodonTransformer" method:
\itemize{
  \item \code{TRUE} (default): Deterministic decoding (most likely tokens)
  \item \code{FALSE}: Probabilistic sampling based on temperature
}}

\item{temperature}{Numeric. Controls randomness in non-deterministic mode 
for "CodonTransformer". Lower values (0.2, default) are conservative; 
higher values (0.8) increase diversity. Must be positive.}

\item{top_p}{Numeric. Nucleus sampling threshold (0-1) for "CodonTransformer". 
Only tokens with cumulative probability up to this value are considered. 
Default: 0.95.}

\item{match_protein}{Logical. For "CodonTransformer", constrains predictions 
to maintain exact amino acid sequence. Recommended for unusual proteins 
or high temperature settings (default: FALSE).}

\item{spliceai}{Logical. Whether to predict splice sites using SpliceAI 
(default: FALSE). Requires appropriate environment setup.}
}
\value{
The return type depends on parameters:
  \itemize{
    \item Single DNAString: When num_sequences=1 and spliceai=FALSE
    \item DNAStringSet: When num_sequences>1 and spliceai=FALSE  
    \item data.table: When spliceai=TRUE (includes sequences and splice predictions)
  }
}
\description{
\code{codon_optimize} redesigns a coding sequence by replacing each codon 
with its optimal synonymous alternative, while maintaining the same amino 
acid sequence. This function supports multiple optimization methods and is 
useful for improving protein expression in heterologous systems.
}
\examples{
cf_all <- count_codons(yeast_cds)
optimal_codons <- est_optimal_codons(cf_all)
seq <- 'ATGCTACGA'
# method "naive":
codon_optimize(seq, optimal_codons)
# method "IDT":
codon_optimize(seq, cf = cf_all, method = "IDT")
codon_optimize(seq, cf = cf_all, method = "IDT", num_sequences = 10)

# # The following examples requires pre-installation of python package SpliceAI or Codon
#  # Transformer. see the codon optimization vignette for further details.
# seq_opt <- codon_optimize(seq, method = "CodonTransformer",
#     organism = "Saccharomyces cerevisiae")
# seqs_opt <- codon_optimize(seq, method = "CodonTransformer",
#     organism = "Saccharomyces cerevisiae", num_sequences = 10,
#     deterministic =FALSE, temperature = 0.4)
# seqs_opt <- codon_optimize(seq, cf = cf_all, method = "IDT",
#     num_sequences = 10, spliceai = TRUE)
# seq_opt <- codon_optimize(seq, method = "CodonTransformer",
#     organism = "Saccharomyces cerevisiae", spliceai = TRUE)
}
\references{
Fallahpour A, Gureghian V, Filion GJ, Lindner AB, Pandi A. CodonTransformer:
a multispecies codon optimizer using context-aware neural networks. Nat Commun. 2025 Apr 3;16(1):3205.

Jaganathan K, Panagiotopoulou S K, McRae J F, et al. Predicting splicing from primary sequence
 with deep learning. Cell, 2019, 176(3): 535-548.e24.
}
