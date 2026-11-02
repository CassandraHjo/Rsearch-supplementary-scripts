library(Rsearch)
fastq_path <- "mock_fastq"

vs_optimize_truncee_rate(file.path(fastq_path, "P15_ZB_R1.fastq.gz"), 
                         file.path(fastq_path, "P15_ZB_R2.fastq.gz"))
