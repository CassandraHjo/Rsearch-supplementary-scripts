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

unoise.tbl <- vs_cluster_unoise(all.tbl,
                                relabel = "ZOTU")

sequence.tbl <- vs_uchime_denovo(unoise.tbl) |> 
  separate_wider_delim(Header, delim = ";", names = c("Header", "size"))

readcount.tbl <- unoise.tbl |> 
  separate_wider_delim(Header, delim = ";", names = c("Header", "size")) |> 
  filter(Header %in% sequence.tbl$Header)
