import os
import subprocess
import argparse

def list_fastq_files(folder_path):
    return [f for f in os.listdir(folder_path) if f.endswith('.fastq.gz')]

def create_pbs_file(fastq_file, job_file_name, input_folder, output_folder, threads, file_type, leading, trailing, slidingwindow, minlen):
    file_id = fastq_file.split('.')[0]
    with open(job_file_name, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("#PBS -l walltime=72:00:00,select=1:ncpus=32:mem=186gb\n")
        f.write(f"#PBS -N Trimmomatic_{file_id}\n")
        f.write("#PBS -A <project_name>\n") #REMOVED FOR PRIVACY
        f.write("#PBS -m abe\n")
        f.write("#PBS -M gdalba@phas.ubc.ca\n")
        f.write("\n")
        f.write("source /project/PROJECT_ENV(REMOVED_FOR_PRIVACY)/gdalba/anaconda3/etc/profile.d/conda.sh\n")
        f.write("source activate /project/PROJECT_ENV(REMOVED_FOR_PRIVACY)/gdalba/anaconda3/envs/trimmomatic\n")
        f.write(f"cd {input_folder}\n")
        

        adapter_file = "TruSeq3-PE.fa" if file_type == "PE" else "TruSeq3-SE.fa"
        trim_options = f"ILLUMINACLIP:/scratch/PROJECT_ENV(REMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RAW_FILES/GENOME/trimmomatic_adapters/{adapter_file}:2:30:10"
        other_options = f"LEADING:{leading} TRAILING:{trailing} SLIDINGWINDOW:{slidingwindow} MINLEN:{minlen}"

        if file_type == "PE":
            paired_file = fastq_file.replace("_1.fastq.gz", "_2.fastq.gz")
            f.write(f"trimmomatic PE -phred33 {fastq_file} {paired_file} {output_folder}/{fastq_file}.paired.trimmed.fastq.gz {output_folder}/{fastq_file}.unpaired.trimmed {output_folder}/{paired_file}.paired.trimmed {output_folder}/{paired_file}.unpaired.trimmed {trim_options} {other_options}\n")

        else:
            f.write(f"trimmomatic SE -phred33 {fastq_file} {output_folder}/{fastq_file}.trimmed.fastq.gz {trim_options} {other_options}\n")

def schedule_job(job_file_name):
    subprocess.run(["qsub", job_file_name])

def main(input_folder, output_folder, threads, file_type, leading, trailing, slidingwindow, minlen):
    fastq_files = list_fastq_files(input_folder)
    
    for fastq_file in fastq_files:
        if file_type == "PE" and "_2.fastq.gz" in fastq_file:
            continue  # Skip the second file in paired-end, as it will be processed with its partner
        
        job_file_name = f"{fastq_file}_trimmomatic.pbs"
        create_pbs_file(fastq_file, job_file_name, input_folder, output_folder, threads, file_type, leading, trailing, slidingwindow, minlen)
        schedule_job(job_file_name)
        print(f"Scheduled Trimmomatic job for {fastq_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate Trimmomatic jobs with PBS")
    parser.add_argument("--input_folder", required=True, help="Path to the folder containing fastq.gz files")
    parser.add_argument("--output_folder", required=True, help="Path to the folder where trimmed files will be saved")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads to use in each PBS job")
    parser.add_argument("--type", choices=["SE", "PE"], required=True, help="Type of fastq files: SE for single-end, PE for paired-end")
    parser.add_argument("--leading", type=int, default=3, help="Quality threshold for LEADING option")
    parser.add_argument("--trailing", type=int, default=3, help="Quality threshold for TRAILING option")
    parser.add_argument("--slidingwindow", default="4:15", help="Parameters for SLIDINGWINDOW option")
    parser.add_argument("--minlen", type=int, default=36, help="Minimum read length to keep")
    args = parser.parse_args()

    main(args.input_folder, args.output_folder, args.threads, args.type, args.leading, args.trailing, args.slidingwindow, args.minlen)