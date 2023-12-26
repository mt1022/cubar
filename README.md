
# cubar

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cubar)](https://CRAN.R-project.org/package=cubar)
[![](https://cranlogs.r-pkg.org/badges/cubar)](https://cran.r-project.org/package=cubar)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10155990.svg)](https://doi.org/10.5281/zenodo.10155990)
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

The latest release of `cubar` can be installed with:

```r
install.packages("cubar")
```

The latest developmental version of `cubar` can be installed with:

```r
devtools::install_github("mt1022/cubar", dependencies = TRUE)
```

### Usage
Documentation can be found within R (by typing `?function_name`). The following tutorials are available from our [website](https://mt1022.github.io/cubar/):

- [Get Started](https://mt1022.github.io/cubar/articles/cubar.html)
- [Non-standard Genetic Code](https://mt1022.github.io/cubar/articles/non_standard_genetic_code.html)
- [Theories behind cubar](https://mt1022.github.io/cubar/articles/theory.html)

### Suggests
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) for sequence input/output and manipulation.
- [Peptides](https://github.com/dosorio/Peptides) for peptide- or protein-related indices.

### Getting help
Please use GitHub [issues](https://github.com/mt1022/cubar/issues) for bug reports, questions, and feature requests.
