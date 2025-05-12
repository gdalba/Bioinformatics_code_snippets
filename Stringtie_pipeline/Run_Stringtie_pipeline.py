import os

paired_end_reads_path = "/scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RAW_FILES/RNA_READS/BROWNE2017/"
paired_end_reads_files = [f for f in os.listdir(paired_end_reads_path) if f.endswith("_R1_paired.fastq.gz")]

for read_file in paired_end_reads_files:
    base_name = read_file.replace("_R1_paired.fastq.gz", "")
    job_name = "assembly_pipeline_" + base_name.replace("wBrowne_Project052012_", "")

    pbs_script = f"""#!/bin/bash

#PBS -l walltime=40:00:00,select=4:ncpus=14:mem=20gb
#PBS -N {job_name}
#PBS -A #REMOVED FOR PRIVACY
#PBS -m abe
#PBS -M gdalba@phas.ubc.ca

source /project/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/anaconda3/etc/profile.d/conda.sh
conda activate hisat2

cd $PBS_O_WORKDIR

hisat2 --dta --repeat --summary-file {base_name}_summary.txt --threads 56 -x /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RAW_FILES/GENOME/Mnemiopsis_leidyi_genome -1 {paired_end_reads_path}/{base_name}_R1_paired.fastq.gz -2 {paired_end_reads_path}/{base_name}_R2_paired.fastq.gz -S {base_name}.sam

conda deactivate
conda activate samtools

samtools view -b -o {base_name}.bam {base_name}.sam
samtools sort -o {base_name}_sorted.bam {base_name}.bam
samtools index {base_name}_sorted.bam

conda deactivate
conda activate stringtie

stringtie -p 56 -l MLE -G {base_name}_sorted.bam -o {base_name}_transcripts.gtf

conda deactivate
conda activate gffread

gffread {base_name}_transcripts.gtf -g /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RAW_FILES/GENOME/Mnemiopsis_leidyi_genome.fa -w {base_name}_transcripts.fa

conda deactivate

conda activate bowtie2

bowtie2-build /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RAW_FILES/GENOME/Mnemiopsis_leidyi_genome.fa Mnemiopsis_leidyi_genome_index

bowtie2 -p 56 -x Mnemiopsis_leidyi_genome_index -1 {paired_end_reads_path}/{base_name}_R1_paired.fastq.gz -2 {paired_end_reads_path}/{base_name}_R2_paired.fastq.gz -S {base_name}_bowtie2.sam

samtools view -b -o {base_name}_bowtie2.bam {base_name}_bowtie2.sam
samtools sort -o {base_name}_bowtie2_sorted.bam {base_name}_bowtie2.bam
samtools index {base_name}_bowtie2_sorted.bam

conda deactivate

conda activate minimap2

minimap2 -ax sr -t 56 /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RAW_FILES/GENOME/Mnemiopsis_leidyi_genome.fa {paired_end_reads_path}/{base_name}_R1_paired.fastq.gz {paired_end_reads_path}/{base_name}_R2_paired.fastq.gz > {base_name}_minimap2.sam

samtools view -b -o {base_name}_minimap2.bam {base_name}_minimap2.sam
samtools sort -o {base_name}_minimap2_sorted.bam {base_name}_minimap2.bam
samtools index {base_name}_minimap2_sorted.bam

conda deactivate


conda activate TransDecoder

TransDecoder.LongOrfs -t {base_name}_transcripts.fa
TransDecoder.Predict -t {base_name}_transcripts.fa --retain_pfam_hits /arc/project/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/pfam/Pfam-A.hmm

conda deactivate
conda activate blast

blastp -query {base_name}_transcripts.fa.transdecoder.pep -db /arc/project/PROJECT_ENVREMOVED_FOR_PRIVACY)/uniref90/uniref90.fasta -num_threads 56 -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -out {base_name}_blastp_uniref90.txt

conda deactivate
conda activate BUSCO546

busco -i {base_name}_transcripts.fa -l /path/to/busco_lineage -o {base_name}_busco -m transcriptome --cpu 56

conda deactivate
conda activate agat

agat_sp_statistics.pl --gff {base_name}_transcripts.gtf -o {base_name}_agat_stats.txt

mv {base_name}_transcripts.gtf /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RESULTS/ANALYSIS/
mv {base_name}_transcripts.fa /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RESULTS/ANALYSIS/
mv {base_name}_blastp_uniref90.txt /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RESULTS/ANALYSIS/
mv {base_name}_busco /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RESULTS/ANALYSIS/
mv {base_name}_agat_stats.txt /scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/RESULTS/ANALYSIS/

rm {base_name}.sam
rm {base_name}.bam
rm {base_name}_sorted.bam
rm {base_name}_sorted.bam.bai
rm -r {base_name}_transcripts.fa.transdecoder_dir
rm {base_name}_transcripts.fa.transdecoder.bed
rm {base_name}_transcripts.fa.transdecoder.gff3
rm {base_name}_transcripts.fa.transdecoder.pep

    """

    with open(f"/scratch/PROJECT_ENVREMOVED_FOR_PRIVACY)/gdalba/PAPER_ANALYSIS/SCRIPTS/PBS/{job_name}.pbs", "w") as pbs_file:
        pbs_file.write(pbs_script)

