R1 <- readFastq(str_c("fastq/", sample.tbl$R1_file[1]))
R2 <- readFastq(str_c("fastq/", sample.tbl$R2_file[1]))

vs_merging_lengths(R1, R2)