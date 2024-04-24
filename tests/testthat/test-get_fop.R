test_that("works for yeast_cds", {
    # determine optimal codons requires sufficient data.
    # we only test the function works.
    set.seed(2024)
    seqs <- yeast_cds[sample.int(length(yeast_cds), 100)]
    expect_visible(get_fop(seqs))
})
