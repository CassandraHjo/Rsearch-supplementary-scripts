for(i in 1:nrow(sample.tbl)){
  R1 <- readFastq(str_c("fastq/", sample.tbl$R1_file[i]))
  R2 <- readFastq(str_c("fastq/", sample.tbl$R2_file[i]))
  vs_fastx_trim_filt(R1, R2, truncee_rate = 0.010) |> 
    vs_fastq_mergepairs(output_format = "fasta") |> 
    vs_fastx_uniques(output_format = "fasta", sample = sample.tbl$sample_id[i]) |> 
    writeFasta(out.file = str_c("tmp/", sample.tbl$sample_id[i], ".fasta"))
}

all.tbl <- lapply(list.files("tmp", full.names = T), readFasta) |> 
  bind_rows()

all_derep.tbl <- vs_fastx_uniques(all.tbl,
                                  output_format = "fasta",
                                  minuniquesize = 2)

all_derep.tbl <- vs_uchime_denovo(all_derep.tbl)

sequence.tbl <- vs_cluster_size(all_derep.tbl, 
                                id = identity,
                                relabel = "OTU")

readcount.tbl <- vs_usearch_global(all.tbl,
                                   database = sequence.tbl,
                                   id = identity,
                                   otutabout = T)
