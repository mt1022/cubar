
# cubar

<!-- badges: start -->
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

## Installation

The latest version of `cubar` can be installed with:

```r
devtools::install_github('mt1022/cubar', dependencies = TRUE)
```

## Example
see `vignettes/Introduction.Rmd`.
