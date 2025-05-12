#!/bin/bash

input_sra_list="SraAccList.txt"

source /project/st-ssplotki-1/gdalba/anaconda3/etc/profile.d/conda.sh
conda activate entrez

echo "SRA accessions to download:"
cat "$input_sra_list"

# Use parallel to download the SRA accessions concurrently
parallel -j 16 'fastq-dump --split-files --gzip {} -O ./' :::: "$input_sra_list"
