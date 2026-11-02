#!/bin/bash
set -euo pipefail

# Settings – lokal Mac
base="dir"
metadata_file_tsv="$base/mock_metadata.txt"
fastq_folder="$base/mock_fastq"
out_folder="$base/vsearch_output"

threads=4

###############################################
### Creating output directories
###

if [ -d "$out_folder" ]
then
    echo "*** Removing existing output folder..."
    rm -rf -- "$out_folder"
fi

echo "*** Creating output directories..."

mkdir -p "$out_folder/filtered"
mkdir -p "$out_folder/merged"
mkdir -p "$out_folder/dereplicated"

############################
### Reading metadata table
###
n_lines=$(awk 'END {print NR}' "$metadata_file_tsv")

###########################
### Looping over samples
###
for i in $(seq 2 "$n_lines")
do
    # Columns: sample_id, run_id, R1_file, R2_file
    read -r SampleID run_id R1_file R2_file <<< "$(sed -n "${i}p" "$metadata_file_tsv")"

    [[ -z "$SampleID" ]] && continue
    R2_file="${R2_file%$'\r'}"

    r1_file="$fastq_folder/$R1_file"
    r2_file="$fastq_folder/$R2_file"

    echo "*** Processing sample $SampleID..."

    # Check input files
    [[ -f "$r1_file" && -f "$r2_file" ]] || {
        echo "Missing FASTQ files for $SampleID"
        exit 1
    }

    ############################
    ### Trimming and filtering
    ###
    echo "*** VSEARCH trimming and filtering..."

    vsearch --fastq_filter "$r1_file" \
        --reverse "$r2_file" \
        --fastq_truncqual 1 \
        --fastq_maxee_rate 0.01 \
        --fastq_minlen 20 \
        --fastq_truncee_rate 0.01 \
        --fastqout "$out_folder/filtered/${SampleID}_R1.fastq" \
        --fastqout_rev "$out_folder/filtered/${SampleID}_R2.fastq"


    ############################
    ### Merging paired reads
    ###
    echo "*** VSEARCH merging R1 and R2..."

    vsearch --fastq_mergepairs "$out_folder/filtered/${SampleID}_R1.fastq" \
        --reverse "$out_folder/filtered/${SampleID}_R2.fastq" \
        --fastqout "$out_folder/merged/${SampleID}_merged.fastq"
    
    ############################
    ### Dereplication
    ###
    vsearch \
    --fastx_uniques "$out_folder/merged/${SampleID}_merged.fastq" \
    --sizeout \
    --sample "$SampleID" \
    --fastaout "$out_folder/dereplicated/${SampleID}_derep.fasta"

done

###############################################
### Writing all reads to the same file
###
echo "*** Combining all samples..."

cat "$out_folder"/dereplicated/*.fasta > "$out_folder/all.fasta"

###############################################
### Global dereplication
###
echo "*** Global dereplication..."

vsearch \
    --fastx_uniques "$out_folder/all.fasta" \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --minuniquesize 2 \
    --fastaout "$out_folder/all_dereplicated.fasta"
    
###############################################
### Global chimera detection
###
echo "*** Detecting chimeras..."

vsearch \
    --uchime_denovo "$out_folder/all_dereplicated.fasta" \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras "$out_folder/all_nonchimeras.fasta"
    
###############################################
### Cluster sequences into OTUs
###
echo "*** Clustering OTUs..."

vsearch \
    --cluster_unoise "$out_folder/all_nonchimeras.fasta" \
    --sizein \
    --sizeout \
    --relabel OTU \
    --threads "$threads" \
    --centroids "$out_folder/OTUs.fasta" \
    --uc "$out_folder/OTUs.uc" \
    --otutabout "$out_folder/read_count.tsv"

echo "*** Finished"
