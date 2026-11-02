library(Rsearch)
fastq_path <- "mock_fastq"

plot_base_quality(file.path(fastq_path, "P15_ZB_R1.fastq.gz"),
                  file.path(fastq_path, "P15_ZB_R2.fastq.gz"),
                  show_overlap_box = FALSE)
