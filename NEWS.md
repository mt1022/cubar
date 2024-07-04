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
