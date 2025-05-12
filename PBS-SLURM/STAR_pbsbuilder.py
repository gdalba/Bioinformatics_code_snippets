import os
import sys
import subprocess
import argparse
from collections import defaultdict

def main(args):
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    fastq_files = os.listdir(args.input)
    base_to_files = defaultdict(list)
    for file_name in fastq_files:
        if file_name.endswith('.fastq.gz'):
            base_name = file_name.rsplit('_', 1)[0] if '_1' in file_name or '_2' in file_name else file_name.rsplit('.', 1)[0]
            base_to_files[base_name].append(os.path.join(args.input, file_name))

    for base_name, file_list in base_to_files.items():
        output_file = os.path.join(args.output, base_name + "_output.txt")
        error_file = os.path.join(args.output, base_name + "_error.txt")
        genome_index_dir = args.genome
        output_prefix = base_name + "_"

        if not os.path.exists(genome_index_dir):
            print("STAR index does not exist. Please generate it before running this script.")
            sys.exit(1)

        if len(file_list) == 1:
            star_cmd = create_star_cmd(args.threads, genome_index_dir, file_list[0], args.output, output_prefix, args.overhang)
        elif len(file_list) == 2:
            star_cmd = create_star_cmd(args.threads, genome_index_dir, file_list[0], args.output, output_prefix, args.overhang, file_list[1])
        else:
            print(f"Unexpected number of files for base name {base_name}: {file_list}")
            continue

        pbs_script = create_pbs_script(args.threads, base_name, output_file, error_file, star_cmd)

        pbs_file_path = os.path.join(args.output, base_name + ".pbs")
        with open(pbs_file_path, "w") as pbs_file:
            pbs_file.write(pbs_script)

        subprocess.run(["qsub", pbs_file_path], check=True)


def create_star_cmd(threads, genome_index_dir, input_file1, output_dir, output_prefix, overhang, input_file2=None):
    if input_file2:
        return f"STAR --runThreadN {threads} --genomeDir {genome_index_dir} --readFilesIn {input_file1} {input_file2} --readFilesCommand zcat --twopassMode Basic --twopass1readsN -1 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {output_dir}/{output_prefix} --sjdbOverhang {overhang} --outSAMtype BAM Unsorted"
    else:
        return f"STAR --runThreadN {threads} --genomeDir {genome_index_dir} --readFilesIn {input_file1} --readFilesCommand zcat --twopassMode Basic --twopass1readsN -1 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {output_dir}/{output_prefix} --sjdbOverhang {overhang} --outSAMtype BAM Unsorted"


def create_pbs_script(threads, analysis_name, output_file, error_file, star_cmd):
    return f"""#!/bin/bash
#PBS -l walltime=36:00:00,select=1:ncpus={threads}:mem=186gb
#PBS -N {analysis_name}
#PBS -A #REMOVED FOR PRIVACY
#PBS -m abe
#PBS -M gdalba@phas.ubc.ca
#PBS -o {output_file}
#PBS -e {error_file}

source /project/PROJECT_ENV(REMOVED_FOR_PRIVACY)/gdalba/anaconda3/etc/profile.d/conda.sh

conda activate STAR

{star_cmd}
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Build PBS scripts for STAR.')
    parser.add_argument('-i', '--input', required=True, help='input folder with .fastq.gz files')
    parser.add_argument('-g', '--genome', required=True, help='directory containing STAR genome index')
    parser.add_argument('-o', '--output', required=True, help='output directory')
    parser.add_argument('-t', '--threads', required=True, type=int, help='number of threads')
    parser.add_argument('-oh', '--overhang', required=True, type=int, help='Junction overhang size')

    args = parser.parse_args()
    main(args)
