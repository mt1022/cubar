library(data.table)
library(rlang)
library(Biostrings)
library(ggplot2)

devtools::load_all("~/codehub/cubar")

cds_seq <- readDNAStringSet('~/Downloads/Homo_sapiens.GRCh38.cds.all.fa.gz')
names(cds_seq) <- sub(' .*', '', names(cds_seq))
#cds_seq <- cds_seq[1:100]

cds_seq <- qc_cds(cds_seq, check_stop = TRUE, rm_start = FALSE)
names(cds_seq) <- sub(' .*', '', names(cds_seq))

set.seed(123)
cds_seq_test <- cds_seq[sample.int(length(cds_seq), 100)]

cds_enc <- get_enc(cds_seq)
fwrite(stack(cds_enc), file = 'sample_enc.tsv', sep = '\t')

get_enc(cds_seq_test)
