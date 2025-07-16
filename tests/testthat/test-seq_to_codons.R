test_that("sequence-to-codons conversion works for a DNAString", {
  # DNAString
  dna_seq <- Biostrings::DNAString("ATGCGT")
  codons <- seq_to_codons(dna_seq)
  expect_equal(codons, c("ATG", "CGT"))
})

test_that("sequence-to-codons conversion works for a single character string", {
  # Single character string
  dna_seq <- "ATGCGT"
  codons <- seq_to_codons(dna_seq)
  expect_equal(codons, c("ATG", "CGT"))
})

test_that("sequence-to-codons conversion works for sequence whose length is not multiple of three", {
  # a sequence with a length of 5
  dna_seq <- Biostrings::DNAString("ATGCG")
  codons <- seq_to_codons(dna_seq)
  expect_equal(codons, c("ATG"))
})

