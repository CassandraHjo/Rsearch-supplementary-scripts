R1 <- readFastq(str_c("fastq/", sample.tbl$R1_file[1]))
R2 <- readFastq(str_c("fastq/", sample.tbl$R2_file[1]))

plot_base_quality(R1, R2, show_overlap_box = F)
