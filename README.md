
# cubar

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

### Overview
cubar is a package for codon usage bias analysis in R. Main features are listed as follows:

- Codon level analyses
    - Support NCBI and custom genetic codes
    - Calculate tRNA weights
    - Calculate relative synonymous codon usage (RSCU)
    - Machine learning-based inference of optimal codons
    - Visualization codon-anticodon pairing relationships
- Gene level analyses
    - Codon frequency matrix
    - Codon Adaptation Index (CAI)
    - Mean Codon Stabilization Coefficients (CSCg)
    - Effective number of codons (ENC)
    - Fraction of optimal codons (Fop)
    - GC content at 4-fold degenerate sites (GC4d)
    - tRNA Adaptation Index (tAI)

### Dependencies
Depends

- `R` (>= 4.1.0)

Imports

- `Biostrings` (>= 2.60.0),
- `IRanges` (>= 2.34.0),
- `data.table` (>= 1.14.0),
- `ggplot2` (>= 3.3.5),
- `rlang` (>= 0.4.11)

### Installation

The latest version of `cubar` can be installed with:

```r
devtools::install_github('mt1022/cubar', dependencies = TRUE)
```

### Suggests
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) for sequence input/output and manipulation.
- [Peptides](https://github.com/dosorio/Peptides) for peptide- or protein-related indices.

### Getting help
Please use GitHub [issues](https://github.com/mt1022/cubar/issues) for bug reports, questions, and feature requests.
