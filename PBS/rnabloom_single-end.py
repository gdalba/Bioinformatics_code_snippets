import os
import sys
import subprocess

print('-----------------------------------')
print('RNABLOOM PBS BUILDER - 06/19/2023')
print('rnabloom_pbsbuilder.py')
print('Made by Gabriel Dall\'Alba')
print('Please provide 3 arguments:')
print('input_folder, out_dir, threads')
print('-----------------------------------')

# Receive sysargs for input folder, output dir, and threads
input_folder = sys.argv[1]
output_dir = sys.argv[2]
threads = sys.argv[3]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

fastq_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.fastq.gz')]

rnabloom_cmd = 'rnabloom -sef {files} -revcomp-right -t {threads} -outdir {outdir}'

pbs_template = """#!/bin/bash
#PBS -l walltime=48:00:00,select=1:ncpus={threads}:mem=186gb
#PBS -N RNABloom_Analysis
#PBS -A PROJECT_ENV(REMOVED_FOR_PRIVACY)
#PBS -m abe
#PBS -M gdalba@phas.ubc.ca
#PBS -o {output_dir}/rnabloom_output.txt
#PBS -e {output_dir}/rnabloom_error.txt

source /project/PROJECT_ENV(REMOVED_FOR_PRIVACY)/gdalba/anaconda3/etc/profile.d/conda.sh

conda activate rnabloom

{rnabloom_cmd}
"""

# Build RNABloom command with all files as input and submit to sockeye
files = ' '.join(fastq_files)
rnabloom_command = rnabloom_cmd.format(files=files, threads=threads, outdir=output_dir)

pbs_script = pbs_template.format(threads=threads, output_dir=output_dir, rnabloom_cmd=rnabloom_command)

pbs_file_path = os.path.join(output_dir, "rnabloom_analysis.pbs")

with open(pbs_file_path, "w") as pbs_file:
    pbs_file.write(pbs_script)

# Submit the PBS script to the job scheduler
subprocess.run(["qsub", pbs_file_path])
