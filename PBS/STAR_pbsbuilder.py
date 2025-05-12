import os
import sys
import subprocess
import argparse
from collections import defaultdict

print('-----------------------------------')
print('CTENOTATION V1 - 07/28/2023')
print('star_pbsbuilder.py')
print('Made by Gabriel Dall\'Alba')
print('Please provide 4 arguments:')
print('-i or --input for input_folder,')
print('-g or --genome for genome_file,')
print('-o or --output for output_dir,')
print('-t or --threads for threads')
print('-oh or --overhang for overhang')
print('-----------------------------------')

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Build PBS scripts with STAR.')
parser.add_argument('-i', '--input', required=True, help='input folder with .fastq.gz files')
parser.add_argument('-g', '--genome', required=True, help='genome fasta file')
parser.add_argument('-o', '--output', required=True, help='output directory')
parser.add_argument('-t', '--threads', required=True, type=int, help='number of threads')
parser.add_argument('-oh', '--overhang', required=True, type=int, help='Overhang parameter for STAR genome Indexing')
args = parser.parse_args()

# Get the arguments
input_folder = args.input
genome_file = args.genome
output_dir = args.output
threads = args.threads
overhang = args.overhang

# Check if output_dir exists, if not create it
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Group fastq files by base name (without _1 or _2 and extension)
fastq_files = os.listdir(input_folder)
base_to_files = defaultdict(list)
for file_name in fastq_files:
    if file_name.endswith('.fastq.gz'):
        base_name = file_name.rsplit('_', 1)[0] if '_1' in file_name or '_2' in file_name else file_name
        base_to_files[base_name].append(os.path.join(input_folder, file_name))

# STAR command for single-end and paired-end
star_cmd_single = "STAR --runThreadN {threads} --genomeDir {genome_index_dir} --readFilesIn {input_files} --readFilesCommand zcat --outFileNamePrefix {output_prefix}"
star_cmd_paired = "STAR --runThreadN {threads} --genomeDir {genome_index_dir} --readFilesIn {input_file1} {input_file2} --readFilesCommand zcat --outFileNamePrefix {output_prefix}"

# PBS script template
pbs_template = """#!/bin/bash
#PBS -l walltime=36:00:00,select=1:ncpus={threads}:mem=186gb
#PBS -N {analysis_name}
#PBS -A (REMOVED_FOR_PRIVACY)
#PBS -m abe
#PBS -M gdalba@phas.ubc.ca
#PBS -o {output_file}
#PBS -e {error_file}

source /project/(REMOVED_FOR_PRIVACY)/gdalba/anaconda3/etc/profile.d/conda.sh

conda activate star

{star_cmd}
"""

# For each base name, create and submit a PBS script
for base_name, file_list in base_to_files.items():
    analysis_name = base_name + "_analysis"
    output_file = os.path.join(output_dir, base_name + "_output.txt")
    error_file = os.path.join(output_dir, base_name + "_error.txt")
    genome_index_dir = os.path.join(output_dir, "star_index")

    if not os.path.exists(genome_index_dir):
        os.makedirs(genome_index_dir)
        subprocess.run(f"STAR --runMode genomeGenerate --genomeDir {genome_index_dir} --genomeFastaFiles {genome_file} --runThreadN {threads} --sjdbOverhang {overhang}", shell=True)

    if len(file_list) == 1:  # Single-end
        star_cmd = star_cmd_single.format(threads=threads, genome_index_dir=genome_index_dir, input_files=file_list[0], output_prefix=base_name + "_")
    elif len(file_list) == 2:  # Paired-end
        star_cmd = star_cmd_paired.format(threads=threads, genome_index_dir=genome_index_dir, input_file1=file_list[0], input_file2=file_list[1], output_prefix=base_name + "_")
    else:
        raise ValueError(f"Unexpected number of files for base name {base_name}: {file_list}")

    pbs_script = pbs_template.format(threads=threads, analysis_name=analysis_name, output_file=output_file, error_file=error_file, star_cmd=star_cmd)

    pbs_file_path = os.path.join(output_dir, base_name + ".pbs")
    
    with open(pbs_file_path, "w") as pbs_file:
        pbs_file.write(pbs_script)

    # Submit the PBS script to the job scheduler
    subprocess.run(["qsub", pbs_file_path])