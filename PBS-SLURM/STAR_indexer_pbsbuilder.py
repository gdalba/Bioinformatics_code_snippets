import os
import argparse
import subprocess

def main(args):
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    star_cmd = create_star_cmd(args.threads, args.genome, args.output)

    pbs_script = create_pbs_script(args.threads, args.output, star_cmd, args.account, 
                                  args.email, args.conda_path, args.conda_env)

    pbs_file_path = os.path.join(args.output, "index_genome.pbs")
    with open(pbs_file_path, "w") as pbs_file:
        pbs_file.write(pbs_script)

    if args.submit:
        subprocess.run(["qsub", pbs_file_path], check=True)
        print(f"Job submitted: {pbs_file_path}")
    else:
        print(f"PBS script created: {pbs_file_path}")
        print("Use 'qsub' command to submit the job")

def create_star_cmd(threads, genome_file, output_dir):
    return f"STAR --runMode genomeGenerate --genomeDir {output_dir} --genomeFastaFiles {genome_file} --runThreadN {threads} --outFileNamePrefix {output_dir}/"

def create_pbs_script(threads, output_dir, star_cmd, account, email, conda_path, conda_env):
    return f"""#!/bin/bash
#PBS -l walltime=48:00:00,select=1:ncpus={threads}:mem=186gb
#PBS -N STAR_indexing
#PBS -A {account}
#PBS -m abe
#PBS -M {email}
#PBS -o {os.path.join(output_dir, "pbs_output.txt")}
#PBS -e {os.path.join(output_dir, "pbs_error.txt")}

source {conda_path}/etc/profile.d/conda.sh

conda activate {conda_env}

{star_cmd}
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Build a PBS script for STAR genome indexing.')
    parser.add_argument('-g', '--genome', required=True, help='genome fasta file')
    parser.add_argument('-o', '--output', required=True, help='output directory for genome index')
    parser.add_argument('-t', '--threads', required=True, type=int, help='number of threads')
    parser.add_argument('-a', '--account', default='your_account_id', help='PBS account allocation')
    parser.add_argument('-e', '--email', default='your.email@example.com', help='email for job notifications')
    parser.add_argument('-c', '--conda-path', default='/path/to/your/anaconda3', 
                      help='path to conda installation')
    parser.add_argument('-v', '--conda-env', default='STAR', help='conda environment with STAR')
    parser.add_argument('-s', '--submit', action='store_true', help='submit the job automatically')
    
    args = parser.parse_args()
    main(args)
