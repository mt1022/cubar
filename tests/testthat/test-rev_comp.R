test_that("rev_comp works for DNAStringSet", {
  # Create a DNAStringSet
  dna <- Biostrings::DNAStringSet(c("ACGT", "TGCA"))
  # Reverse complement
  dna_rev <- rev_comp(dna)
  # Check the result
  expect_equal(as.character(dna_rev), as.character(dna))
})

test_that("rev_comp works for character vector of sequences", {
  # Create a DNAStringSet
  dna <- c("ACGT", "TGCA")
  # Reverse complement
  dna_rev <- rev_comp(dna)
  # Check the result
  expect_equal(as.character(dna_rev), dna)
})
