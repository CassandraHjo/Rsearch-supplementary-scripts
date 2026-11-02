# Load libraries
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(microseq)))
suppressWarnings(suppressMessages(library(dada2)))

# Setup paths
fastq.folder <- "mock_fastq"
metadata.file <- "mock_metadata.txt"
tmp.folder <- "tmp_dada2"
dada2_dir <- "dada2_results"

# Load metadata
metadata.tbl <- read_tsv(metadata.file, show_col_types = FALSE)

# Extract sample names
sample_names <- metadata.tbl$sample_id

# Threads
threads <- 4

# Trimming 

## Prepare outputnames for filtered reads
filt_R1_files <- file.path(tmp.folder, paste0(metadata.tbl$sample_id, "_R1_filt.fastq"))
filt_R2_files <- file.path(tmp.folder, paste0(metadata.tbl$sample_id, "_R2_filt.fastq"))

## Defining arguments
truncLen <- c(250, 220)

# Running filtering
filt.tbl <- filterAndTrim(fwd = file.path(fastq.folder, metadata.tbl$R1_file), 
                          filt = filt_R1_files, 
                          rev = file.path(fastq.folder, metadata.tbl$R2_file), 
                          filt.rev = filt_R2_files,
                          truncLen = truncLen, 
                          compress = FALSE,
                          multithread = threads
)

# Learn errors

err.obj_R1 <- learnErrors(fls = filt_R1_files, 
                          multithread = threads,
                          verbose= TRUE)

err.obj_R2 <- learnErrors(fls = filt_R2_files, 
                          multithread = threads,
                          verbose = TRUE)

# DADA2 main inference

dada.obj_R1 <- dada(filt_R1_files, 
                    err = err.obj_R1, 
                    multithread = threads)

dada.obj_R2 <- dada(filt_R2_files, 
                    err = err.obj_R2, 
                    multithread = threads)

# Merge sequences

merged.obj <- mergePairs(dadaF = dada.obj_R1,
                         derepF = filt_R1_files,
                         dadaR = dada.obj_R2,
                         derepR = filt_R2_files)

# Make read count matrix

readcount.mat <- makeSequenceTable(merged.obj)

# Remove chimeras

readcount.mat <- removeBimeraDenovo(readcount.mat, 
                                    method = "consensus", 
                                    multithread = threads)
