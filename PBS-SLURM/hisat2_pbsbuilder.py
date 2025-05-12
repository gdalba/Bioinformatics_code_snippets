import os
import sys
import subprocess
from collections import defaultdict

def print_header():
    print('-----------------------------------')
    print('HISAT2 PBS Script Generator')
    print('-----------------------------------')
    print('Please provide 4 arguments:')
    print('input_folder, genome_file, out_dir, threads')
    print('-----------------------------------')

print_header()

if len(sys.argv) != 5:
    print(f"Error: Expected 4 arguments, got {len(sys.argv)-1}")
    print("Usage: python hisat2_pbsbuilder.py input_folder genome_file output_dir threads")
    sys.exit(1)

input_folder = sys.argv[1]
genome_file = sys.argv[2]
output_dir = sys.argv[3]
threads = sys.argv[4]

if not os.path.exists(input_folder):
    print(f"Error: Input folder '{input_folder}' does not exist.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")

fastq_files = os.listdir(input_folder)
base_to_files = defaultdict(list)
for file_name in fastq_files:
    if file_name.endswith('.fastq.gz'):
        base_name = file_name.rsplit('_', 1)[0] if '_1' in file_name or '_2' in file_name else file_name
        base_to_files[base_name].append(os.path.join(input_folder, file_name))

hisat2_cmd_single = "hisat2 -x {genome_index} -U {input_files} --dta --sensitive --qc-filter --summary-file {summary_file} -S {sam_file} -p {threads}"
hisat2_cmd_paired = "hisat2 -x {genome_index} -1 {input_file1} -2 {input_file2} --dta --sensitive --qc-filter --summary-file {summary_file} -S {sam_file} -p {threads}"

pbs_template = """#!/bin/bash
#PBS -l walltime=36:00:00,select=1:ncpus={threads}:mem=32gb
#PBS -N {analysis_name}
#PBS -m abe
#PBS -o {output_file}
#PBS -e {error_file}

{hisat2_cmd}

samtools view -@ {threads} -bS {sam_file} > {bam_file}

rm {sam_file}

samtools sort -@ {threads} {bam_file} -o {sorted_bam_file}

rm {bam_file}
"""

print(f"Found {len(base_to_files)} sample(s) to process.")

for base_name, file_list in base_to_files.items():
    analysis_name = base_name + "_analysis"
    output_file = os.path.join(output_dir, base_name + "_output.txt")
    error_file = os.path.join(output_dir, base_name + "_error.txt")
    summary_file = os.path.join(output_dir, base_name + "_summary.txt")
    sam_file = os.path.join(output_dir, base_name + ".sam")
    bam_file = os.path.join(output_dir, base_name + ".bam")
    sorted_bam_file = os.path.join(output_dir, base_name + "_sorted.bam")
    genome_index = genome_file

    if len(file_list) == 1:
        hisat2_cmd = hisat2_cmd_single.format(genome_index=genome_index, input_files=file_list[0], 
                                             summary_file=summary_file, sam_file=sam_file, threads=threads)
    elif len(file_list) == 2:
        hisat2_cmd = hisat2_cmd_paired.format(genome_index=genome_index, input_file1=file_list[0], 
                                             input_file2=file_list[1], summary_file=summary_file, 
                                             sam_file=sam_file, threads=threads)
    else:
        print(f"Warning: Unexpected number of files for base name {base_name}: {file_list}")
        continue

    pbs_script = pbs_template.format(threads=threads, analysis_name=analysis_name, output_file=output_file, error_file=error_file, 
                                   hisat2_cmd=hisat2_cmd, sam_file=sam_file, bam_file=bam_file, 
                                   sorted_bam_file=sorted_bam_file)

    pbs_file_path = os.path.join(output_dir, base_name + ".pbs")
    
    with open(pbs_file_path, "w") as pbs_file:
        pbs_file.write(pbs_script)
    
    print(f"Created PBS script: {pbs_file_path}")

    try:
        subprocess.run(["qsub", pbs_file_path], check=True)
        print(f"Submitted job for {base_name}")
    except (subprocess.SubprocessError, FileNotFoundError) as e:
        print(f"Failed to submit job for {base_name}: {e}")
        print("You may need to submit the job manually.")

print("All PBS scripts have been created and processed.")
