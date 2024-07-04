test_that("differential testing are correct", {
    x <- paste0(rep(c('TTT', 'TTC'), c(10, 20)))
    y <- paste0(rep(c('TTT', 'TTC'), c(20, 10)))

    cudiff <- codon_diff(x, y)
    cudiff <- cudiff[codon %in% c('TTT', 'TTC')]

    f <- fisher.test(rbind(c(10, 20), c(20, 10)))

    expect_true(all(cudiff$global_p == f$p.value))
    expect_true(all(cudiff$fam_p == f$p.value))
    expect_true(all(cudiff$subfam_p == f$p.value))

})
