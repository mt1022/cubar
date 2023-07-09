
# cubar

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Codon usage bias analysis using R!

Features

- Codon level analyses
    - Support NCBI and custom genetic codes
    - Calculate tRNA weights
    - Calculate relative synonymous codon usage (RSCU)
    - Machine learning analyses of preferred codons
    - Show possible codon-anticodon pairings
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

Full Documentation can be found at cubar [website](https://mt1022.github.io/cubar/).


### Suggests
- [Peptides](https://github.com/dosorio/Peptides) for peptide- or protein-related indices.

