
# cubar

> **Comprehensive Codon Usage Bias Analysis in R**

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cubar)](https://CRAN.R-project.org/package=cubar)
[![](https://cranlogs.r-pkg.org/badges/cubar)](https://cran.r-project.org/package=cubar)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10155990.svg)](https://doi.org/10.5281/zenodo.10155990)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

<!-- badges: end -->

## Table of Contents

- [Overview](#overview)
- [Features](#features)
  - [🧬 Codon-Level Analysis](#-codon-level-analysis)
  - [📊 Gene-Level Metrics](#-gene-level-metrics)
  - [🛠️ Utilities \& Tools](#️-utilities--tools)
- [Why Choose cubar?](#why-choose-cubar)
- [Installation](#installation)
  - [Stable Release (Recommended)](#stable-release-recommended)
  - [Development Version](#development-version)
  - [Dependencies](#dependencies)
- [Documentation \& Tutorials](#documentation--tutorials)
  - [🎯 Getting Started](#-getting-started)
  - [📚 Advanced Topics](#-advanced-topics)
- [Example Workflow](#example-workflow)
- [🆘 Getting Help](#-getting-help)
- [Related Packages](#related-packages)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [Citation](#citation)

## Overview

Codon usage bias refers to the non-uniform usage of synonymous codons (codons that encode the same amino acid) across different organisms, genes, and functional categories. **cubar** is a comprehensive R package for analyzing codon usage bias in coding sequences. It provides a unified framework for calculating established codon usage metrics, conducting sliding-window analyses or differential usage analyses, and optimizing sequences for heterologous expression.


## Features

### 🧬 Codon-Level Analysis
- **RSCU calculation**: Relative synonymous codon usage analysis
- **Amino acid usage**: Frequency of each amino acid in sequences
- **Codon weights**: Calculate weights based on gene expression, tRNA availability, and mRNA stability
- **Optimal codon inference**: Machine learning-based identification of optimal codons
- **Codon-anticodon visualization**: Visualization of codon-tRNA pairing relationships

### 📊 Gene-Level Metrics  
- **Codon frequency tabulation**: Count codon occurrences across sequences
- **CAI (Codon Adaptation Index)**: Measure similarity to highly expressed genes 
- **ENC (Effective Number of Codons)**: Assess codon usage bias strength
- **Fop (Fraction of Optimal codons)**: Calculate proportion of optimal codons
- **tAI (tRNA Adaptation Index)**: Match codon usage to tRNA availability
- **CSCg (Codon Stabilization Coefficients)**: Quantify mRNA stability effects 
- **Dp (Deviation from Proportionality)**: Analyze virus-host codon usage relationships
- **GC content metrics**: Overall GC, GC3s (3rd codon positions), GC4d (4-fold degenerate sites)

### 🛠️ Utilities & Tools
- **Sliding window analysis**: Positional codon usage patterns within genes
- **Sequence optimization**: Redesign sequences for optimal expression
- **Differential codon usage**: Statistical comparison between sequence sets
- **Quality control**: Comprehensive CDS validation and preprocessing


## Why Choose cubar?

- **🚀 High Performance**: Process large datasets (>100,000 sequences) efficiently using optimized `Biostrings` and `data.table` backends
- **🧬 Flexible Genetic Codes**: Support for all NCBI genetic codes plus custom genetic code tables
- **🔗 R Ecosystem Integration**: Seamlessly integrate with other bioinformatics and data analysis packages
- **📚 Comprehensive Documentation**: Extensive tutorials, examples, and theoretical background
- **🔬 Research Ready**: Implements established metrics with proper citations and validation


## Installation

### Stable Release (Recommended)

Install the latest stable version from CRAN:

```r
install.packages("cubar")
```

### Development Version

Install the latest development version from GitHub:

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install cubar from GitHub
devtools::install_github("mt1022/cubar", dependencies = TRUE)
```

### Dependencies

**System Requirements:**
- R (≥ 4.1.0)

**Required Packages:**
- `Biostrings` (≥ 2.60.0) - Bioconductor package for sequence manipulation
- `IRanges` (≥ 2.34.0) - Bioconductor infrastructure for range operations  
- `data.table` (≥ 1.14.0) - High-performance data manipulation
- `ggplot2` (≥ 3.3.5) - Data visualization
- `rlang` (≥ 0.4.11) - Language tools

**Note:** Bioconductor packages will be installed automatically, but you may need to update your R installation if you encounter compatibility issues.

## Documentation & Tutorials
📖 **Complete documentation** is available within R (`?function_name`) and on our [**package website**](https://mt1022.github.io/cubar/).

### 🎯 Getting Started
- [**Introduction to cubar**](https://mt1022.github.io/cubar/articles/cubar.html) - Basic usage and core functionality
- [**Non-standard Genetic Codes**](https://mt1022.github.io/cubar/articles/non_standard_genetic_code.html) - Working with alternative genetic codes
- [**Codon Optimization**](https://mt1022.github.io/cubar/articles/codon_optimization.html) - Sequence optimization strategies

### 📚 Advanced Topics  
- [**Mathematical Foundations**](https://mt1022.github.io/cubar/articles/theory.html) - Detailed theory behind the metrics
- [**Function Reference**](https://mt1022.github.io/cubar/reference/) - Complete function documentation

## Example Workflow

Here's a typical analysis workflow demonstrating key functionality:

```r
library(cubar)
library(ggplot2)

# 1. Load and quality-check sequences
data(yeast_cds)
clean_cds <- check_cds(yeast_cds)

# 2. Calculate codon frequencies
codon_freq <- count_codons(clean_cds)

# 3. Calculate multiple metrics
enc <- get_enc(codon_freq)           # Effective number of codons
gc3s <- get_gc3s(codon_freq)         # GC content at 3rd positions

# 4. Analyze highly expressed genes
data(yeast_exp)
yeast_exp <- yeast_exp[yeast_exp$gene_id %in% rownames(codon_freq), ]
high_expr <- head(yeast_exp[order(-yeast_exp$fpkm), ], 500)
rscu_high <- est_rscu(codon_freq[high_expr$gene_id, ])
cai <- get_cai(codon_freq, rscu_high)

# 5. Visualize results
df <- data.frame(ENC = enc, CAI = cai, GC3s = gc3s)
ggplot(df, aes(color = GC3s, x = ENC, y = CAI)) + 
  geom_point(alpha = 0.6) + 
  scale_color_viridis_c() +
  labs(title = "Codon Usage Bias Relationships",
       x = "Effective Number of Codons", y = "Codon Adaptation Index")
```

## 🆘 Getting Help

- **📋 GitHub Issues**: [Report bugs, request features, or ask questions](https://github.com/mt1022/cubar/issues)
- **📖 Documentation**: Check function help (`?function_name`) and [online docs](https://mt1022.github.io/cubar/)


## Related Packages
For complementary analysis, consider these R packages:

- **[Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)** - Sequence input/output and manipulation
- **[Peptides](https://github.com/dosorio/Peptides)** - Peptide and protein property calculations  


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **GitHub Copilot** was used to suggest code snippets during development
- **[GitHub Education](https://education.github.com/)** for providing free access to development tools
- The R and Bioconductor communities for excellent foundational packages
- Contributors and users who have provided feedback and improvements

## Citation
If you use cubar in your research, please cite:

> Mengyue Liu, Bu Zi, Hebin Zhang, Hong Zhang, cubar: a versatile package for codon usage bias analysis in R, Genetics, 2025, iyaf191, https://doi.org/10.1093/genetics/iyaf191

Please also cite the original studies associated with each codon usage metric or third-party software. You can find the relevant references in the documentation of the corresponding functions (for example, type `?cubar::get_enc` in the R console and check the "References" section in the help page).

---

<div align="center">

**[📚 Documentation](https://mt1022.github.io/cubar/) • [🐛 Report Bug](https://github.com/mt1022/cubar/issues) • [💡 Request Feature](https://github.com/mt1022/cubar/issues)**

</div>
