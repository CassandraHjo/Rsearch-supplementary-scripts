# Load libraries
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(microseq)))
suppressWarnings(suppressMessages(library(Rsearch)))

# Setup paths
fastq.folder <- "mock_fastq"
metadata.file <- "mock_metadata.txt"

# Load metadata
metadata.tbl <- read_tsv(metadata.file, show_col_types = FALSE)

# Initiate empty table
all.tbl <- tibble()

# Threads
threads <- 4

# Run preprocessing for each sample
for(i in 1:nrow(metadata.tbl)){
  
  R1.tbl <- readFastq(file.path(fastq.folder, metadata.tbl$R1_file[i]))
  R2.tbl <- readFastq(file.path(fastq.folder, metadata.tbl$R2_file[i]))
  
  # Trim and filter
  R1_filt.tbl <- vs_fastx_trim_filt(file.path(fastq.folder, metadata.tbl$R1_file[i]), 
                                    file.path(fastq.folder, metadata.tbl$R2_file[i]),
                                    minlen = 20, 
                                    truncee_rate = 0.01, 
                                    maxee_rate = 0.01
  )
  R2_filt.tbl <- attr(R1_filt.tbl, "reverse")
  
  # Merge reads
  merged.tbl <- vs_fastq_mergepairs(R1_filt.tbl, R2_filt.tbl, threads = threads)
  
  # Dereplicate
  derep.tbl <- vs_fastx_uniques(fastx_input = merged.tbl,
                                sample = metadata.tbl$sample_id[i])
  
  # Add results to table
  all.tbl <- bind_rows(all.tbl, derep.tbl)
}

# Dereplicate all reads together
all_derep.tbl <- vs_fastx_uniques(all.tbl, minuniquesize = 2)

# Global chimera detection
all_derep_wo_chimera.tbl <- vs_uchime_denovo(all_derep.tbl)

# Clustering
sequence.tbl <- vs_cluster_unoise(all_derep_wo_chimera.tbl,
                                  relabel = "OTU",
                                  threads = threads)
