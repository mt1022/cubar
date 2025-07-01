test_that("works for internal yeast_trna_gcn", {
  w <- est_trna_weight(yeast_trna_gcn)
  expect_equal(dim(w), c(61, 9))
})

test_that('works for a vector of only two anticodons', {
  w <- est_trna_weight(c('Phe-AAA' = 1, 'Phe-GAA' = 1))
  expect_equal(w$W[1:2], c(1.2139, 1.5341))
})

test_that("tai weight values calculated by cubar and the tAI package are equal in yeast.", {
  cubar_w <- est_trna_weight(yeast_trna_gcn, s = list(IC = 0.28, IA = 0.9999, GU = 0.41, UG = 0.68, IU = 0))
  # tai weight values calculated by the tAI package can be calculated using the following code:
  # codon_table = get_codon_table()
  # codon_table[, anticodon := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(codon_table$codon)))]
  # codon_table[, ac_level := yeast_trna_gcn[paste(amino_acid, anticodon, sep = '-')]]
  # codon_table$ac_level[is.na(codon_table$ac_level)] <- 0
  # tai_w <- tAI::get.ws(tRNA=codon_table$ac_level, sking=0)
  tai_w <- c(0.363300492610838, 0.61576354679803, 0.431034482758621, 0.753694581280788,
             0.677339901477833, 0.487684729064039, 0.184796798029557, 0.120689655172414,
             0.29064039408867, 0.492610837438424, 0.145320197044335, 0.246305418719212,
             0.369458128078818, 0.0363300492610838, 0.061576354679803, 0.184729064039409,
             0.0591133004926108, 0.123152709359606, 0.0886699507389163, 0.615775862068966,
             0.197044334975369, 0.254310344827586, 0.431034482758621, 0.554187192118227,
             0.238916256157635, 0.369458128078818, 0.266009852216749, 3.69458128078777e-05,
             0.061576354679803, 0.800492610837439, 0.576354679802956, 0.12323275862069,
             0.677339901477833, 0.487684729064039, 0.24637315270936, 0.140394088669951,
             0.363300492610838, 0.61576354679803, 0.431034482758621, 1, 0.0726600985221675,
             0.123152709359606, 0.677339901477833, 0.278325123152709, 0.862068965517241,
             0.620689655172414, 0.123238916256158, 0.16256157635468, 0.677339901477833,
             0.487684729064039, 0.307949507389163, 0.0985221674876847, 0.58128078817734,
             0.985221674876847, 0.862068965517241, 0.399014778325123, 0.58128078817734,
             0.985221674876847, 0.184729064039409, 0.182266009852217)
  expect_equal(cubar_w$w[-33], tai_w)
})
