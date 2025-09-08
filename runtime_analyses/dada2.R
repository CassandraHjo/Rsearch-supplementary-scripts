filt_R1_file <- file.path("tmp", str_c("filtered_", sample.tbl$R1_file))
filt_R2_file <- file.path("tmp", str_c("filtered_", sample.tbl$R2_file))
names(filt_R1_file) <- sample.tbl$sample_id
names(filt_R2_file) <- sample.tbl$sample_id
flt.tbl <- filterAndTrim(file.path("fastq", sample.tbl$R1_file), filt_R1_file,
                         file.path("fastq", sample.tbl$R2_file), filt_R2_file,
                         truncLen = c(230,200),
                         maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                         compress = F, multithread = F)

err_R1 <- learnErrors(filt_R1_file, verbose = T, multithread = F)
err_R2 <- learnErrors(filt_R2_file, verbose = T, multithread = F)
dada_R1 <- dada(filt_R1_file, err_R1, multithread = F)
dada_R2 <- dada(filt_R2_file, err_R2, multithread = F)

merged.lst <- mergePairs(dada_R1, filt_R1_file,
                         dada_R2, filt_R2_file,
                         verbose = F)
readcount.mat <- makeSequenceTable(merged.lst)
readcount.mat <- removeBimeraDenovo(readcount.mat,
                                    method = "consensus",
                                    multithread = F,
                                    verbose = T)
sequence.tbl <- tibble(Header = str_c("ASV_", 1:ncol(readcount.mat)),
                       Sequence = colnames(readcount.mat),
                       otu_id = Header)
colnames(readcount.mat) <- sequence.tbl$Header
readcount.tbl <- t(readcount.mat) |>
  as_tibble(rownames = "otu_id")
