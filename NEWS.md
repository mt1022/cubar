# cubar 1.0.0

* `est_optimal_codons` and `get_fop` now work on codon frequency matrix like
  other cubar functions.
  
* codon optimization can be done at both family(amino acid) or subfamily level
  now and optimal codons can be estimated for each level using either codon
  bias or gene expression levels (Thanks @maltesemike for valuable suggestions
  and feedback). The false discovery rate is controlled by the `fdr` argument.
  
* There were two RSCU columns (`RSCU` and `rscu`) in the output of
  `est_optimal_codons` and `get_fop`. Now only `rscu` is kept and represents
  the RSCU values.
  
* New functions to perform sliding window analysis on codon usage: `slide`,
  `slide_codon`, `slide_apply` and `slide_plot`.

* New function to calculate the deviation from proportionality (Dp) of host
  tRNA availability: `get_dp`.

# cubar 0.6.0

* Add util functions (`codon_optimize` & `codon_diff`)
* Fix url failure of gtRNAdb, which caused remove of cubar from cran :(

# cubar 0.5.1

* fix a bug in `get_cscg` that caused an error when the input codon frequency
  matrix has a single row.
* finish unit tests for all functions and internal data.

# cubar 0.5.0

* fixed a bug in `est_trna_weight`. Now zero w values were replaced with
  geometric mean (rather than the arithmetic mean) of non-zero w values.
* fixed unexpected warnings in `est_optimal_codons`.
* fixed bugs that update input codon table due to `data.table` reference
  semantics.
* added a new vignette explaining the mathematical details of implementation.

# cubar 0.4.2

* adjust formatting
* fix a typo in `get_enc` code

# cubar 0.4.1

* New vignette for mitochondrial codon usage analysis.
* Fix a bug when in `get_enc` for non-standard genetic code.

# cubar 0.4.0

* Released to CRAN.

# cubar 0.3.2

* Initial CRAN submission.
