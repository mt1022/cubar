
# cubar

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cubar)](https://CRAN.R-project.org/package=cubar)
[![](https://cranlogs.r-pkg.org/badges/cubar)](https://cran.r-project.org/package=cubar)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10155990.svg)](https://doi.org/10.5281/zenodo.10155990)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

### Overview
cubar is a package for codon usage bias analysis in R. Main features are as follows:

- Codon level analyses
    - Calculate tRNA weights;
    - Calculate relative synonymous codon usage (RSCU);
    - Machine learning-based inference of optimal codons;
    - Visualization codon-anticodon pairing relationships;
- Gene level analyses
    - Tabulate codon frequency of each coding sequence;
    - Measure codon usage similarity to highly expressed genes with Codon Adaptation Index (CAI);
    - Quantify the influnce of codon usage on mRNA stability with Mean Codon Stabilization Coefficients (CSCg);
    - Measure codon usage bias with the nonparametric index Effective number of codons (ENC);
    - Measure the fraction of pre-determined optimal codons (Fop) in each sequence;
    - Overall GC content (GC) or that of 3rd synonymous positions (GC3s) or 4-fold degenerate sites (GC4d);
    - Quantify whether codon usage matches tRNA availability using tRNA Adaptation Index (tAI);
- Utilities
    - Sliding window analysis of codon usage within a coding sequence;
    - Optimize codon usage based on optimal codons for heterologous expression;
    - Test differential usage of codons between two sets of sequences;

Main advantages of `cubar` are as follows:
- Process large datasets (>10,0000 sequences) efficiently using the `Biostrings` and `data.table` backends;
- Support genetic codes cataloged by [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) as well as custom ones;
- Integrate with other data analysis or bioinformatic packages in the R ecosystem;

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

- [Get Started](https://mt1022.github.io/cubar/articles/cubar.html): A brief introduction demonstrating the basic usage of `cubar`;
- [Non-standard Genetic Code](https://mt1022.github.io/cubar/articles/non_standard_genetic_code.html): How to use `cubar` with non-standard genetic codes;
- [Theories behind cubar](https://mt1022.github.io/cubar/articles/theory.html): The mathematical details behind the core functions in `cubar`;

### Getting help
Please use GitHub [issues](https://github.com/mt1022/cubar/issues) for bug reports, questions, and feature requests.

### Suggests
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) for sequence input/output and manipulation;
- [Peptides](https://github.com/dosorio/Peptides) for peptide- or protein-related indices;

### Acknowledgements
GitHub Copilot was used to suggest code snippets in the development of this package. Thanks the [GitHub Education](https://education.github.com/) teacher program for providing free access to GitHub Copilot.
