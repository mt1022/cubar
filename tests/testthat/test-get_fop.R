test_that("works for yeast_cds", {
    # determine optimal codons requires sufficient data.
    # we only test the function works.
    set.seed(2024)
    seqs_sample <- yeast_cds[sample.int(length(yeast_cds), 100)]
    cf_sample <- count_codons(seqs_sample)
    expect_visible(get_fop(cf_sample))
})
