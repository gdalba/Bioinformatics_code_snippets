#!/usr/bin/env python3
"""

Author: Gabriel Dall\'Alba
Date: May 10 2025
Version: 1.0

STAR RNA-Seq Alignment Pipeline

This script automates STAR alignment workflow for RNA-seq data. It handles:
- Setting up conda environment with necessary tools
- Building STAR genome index if needed
- Running alignment for single or paired-end data
- Processing all FASTQ files in a directory

Requirements:
- Conda (Miniconda or Anaconda)
- Python 3.6+

Setup Instructions:
1. Make sure conda is installed
2. Run this script with required parameters
3. The script will create a conda environment with STAR and samtools

Usage:
    python run_STAR.py --genome /path/to/genome.fa --fastq_dir /path/to/fastq_files
"""

import os
import sys
import argparse
import subprocess
import logging
import math
import glob
from typing import Dict, List, Tuple

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('star_alignment.log')
    ]
)
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run STAR alignment on fastq files.')
    parser.add_argument('-g', '--genome', required=True, help='Path to reference genome FASTA file')
    parser.add_argument('-f', '--fastq_dir', required=True, help='Directory containing FASTQ files')
    parser.add_argument('-o', '--output_dir', default='star_output', help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads for STAR')
    parser.add_argument('--star_index_dir', help='STAR index directory (optional)')
    parser.add_argument('--conda_env', default='star_env', help='Conda environment name')
    parser.add_argument('--force_index', action='store_true', help='Force rebuild STAR index')
    return parser.parse_args()

def validate_inputs(genome_path: str, fastq_dir: str) -> bool:
    """Validate input files and directories exist."""
    if not os.path.isfile(genome_path):
        logger.error(f"Genome file not found: {genome_path}")
        return False
        
    if not os.path.isdir(fastq_dir):
        logger.error(f"FASTQ directory not found: {fastq_dir}")
        return False
        
    fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq*")) + \
                  glob.glob(os.path.join(fastq_dir, "*.fq*"))
                 
    if not fastq_files:
        logger.error(f"No FASTQ files found in directory: {fastq_dir}")
        return False
        
    return True

def setup_conda_environment(env_name: str) -> bool:
    """Set up conda environment with required packages."""
    try:
        # Check if conda is installed
        subprocess.run(['conda', '--version'], check=True, stdout=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("Conda not found. Please install Miniconda or Anaconda first.")
        return False
        
    # Check if environment exists
    result = subprocess.run(
        f"conda env list | grep '{env_name}' | wc -l",
        shell=True, capture_output=True, text=True
    )
    
    if int(result.stdout.strip()) == 0:
        logger.info(f"Creating conda environment: {env_name}")
        try:
            # Create conda environment with STAR and samtools
            subprocess.run(
                f"conda create -y -n {env_name} -c bioconda star samtools",
                shell=True, check=True
            )
            logger.info(f"Successfully created conda environment: {env_name}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create conda environment: {e}")
            return False
    else:
        logger.info(f"Using existing conda environment: {env_name}")
    
    return True

def get_strandedness() -> str:
    """Ask user if data is single-end or paired-end."""
    while True:
        response = input("Is your data single-end or paired-end? (single/paired): ").lower()
        if response in ['single', 's']:
            return "single"
        elif response in ['paired', 'p']:
            return "paired"
        else:
            print("Invalid input. Please enter 'single' or 'paired'.")

def build_star_index(genome_path: str, index_dir: str, threads: int, conda_env: str) -> bool:
    """Build STAR index from reference genome."""
    os.makedirs(index_dir, exist_ok=True)
    
    # Calculate appropriate genomeSAindexNbases based on genome size
    genome_size = os.path.getsize(genome_path)
    genomeSAindexNbases = min(14, int(math.log2(genome_size)/2) - 1)
    
    logger.info(f"Building STAR index in {index_dir}")
    try:
        cmd = [
            "conda", "run", "-n", conda_env, 
            "STAR", "--runMode", "genomeGenerate",
            "--genomeDir", index_dir,
            "--genomeFastaFiles", genome_path,
            "--runThreadN", str(threads),
            "--genomeSAindexNbases", str(genomeSAindexNbases)
        ]
        
        logger.info(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        logger.info("STAR index built successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to build STAR index: {e}")
        return False

def find_fastq_files(fastq_dir: str) -> List[str]:
    """Find all FASTQ files in the directory."""
    fastq_files = []
    for ext in ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]:
        fastq_files.extend(glob.glob(os.path.join(fastq_dir, ext)))
    return fastq_files

def group_paired_fastq_files(fastq_files: List[str]) -> Dict[str, Tuple[str, str]]:
    """Group paired-end FASTQ files based on naming patterns."""
    paired_files = {}
    
    # Track potential pairs
    r1_files = {}
    r2_files = {}
    
    for file_path in fastq_files:
        basename = os.path.basename(file_path)
        
        # Handle common naming patterns for paired files
        sample_id = None
        
        # Pattern: _R1/_R2
        if "_R1" in basename:
            sample_id = basename.split("_R1")[0]
            r1_files[sample_id] = file_path
        elif "_R2" in basename:
            sample_id = basename.split("_R2")[0]
            r2_files[sample_id] = file_path
            
        # Pattern: _1/_2
        elif "_1." in basename:
            sample_id = basename.split("_1.")[0]
            r1_files[sample_id] = file_path
        elif "_2." in basename:
            sample_id = basename.split("_2.")[0]
            r2_files[sample_id] = file_path
    
    # Match R1 and R2 files
    for sample_id in r1_files:
        if sample_id in r2_files:
            paired_files[sample_id] = (r1_files[sample_id], r2_files[sample_id])
    
    return paired_files

def run_star_alignment(index_dir: str, fastq_files: List[str], output_dir: str, 
                      strandedness: str, threads: int, conda_env: str) -> bool:
    """Run STAR alignment on fastq files."""
    os.makedirs(output_dir, exist_ok=True)
    
    if strandedness == "single":
        # Process single-end fastq files
        for fastq_file in fastq_files:
            sample_name = os.path.basename(fastq_file).split(".")[0]
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)
            
            logger.info(f"Running STAR alignment for sample: {sample_name}")
            try:
                cmd = [
                    "conda", "run", "-n", conda_env,
                    "STAR", "--genomeDir", index_dir,
                    "--readFilesIn", fastq_file,
                    "--outFileNamePrefix", f"{sample_output_dir}/",
                    "--runThreadN", str(threads),
                    "--outSAMtype", "BAM", "SortedByCoordinate",
                    "--outSAMattributes", "Standard"
                ]
                
                # Handle gzipped files
                if fastq_file.endswith(".gz"):
                    cmd.extend(["--readFilesCommand", "zcat"])
                
                subprocess.run(cmd, check=True)
                logger.info(f"Alignment completed for {sample_name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to run alignment for {sample_name}: {e}")
                return False
    else:
        # Process paired-end fastq files
        paired_files = group_paired_fastq_files(fastq_files)
        
        if not paired_files:
            logger.error("No valid paired-end FASTQ file pairs found.")
            return False
        
        for sample_id, (r1_file, r2_file) in paired_files.items():
            sample_output_dir = os.path.join(output_dir, sample_id)
            os.makedirs(sample_output_dir, exist_ok=True)
            
            logger.info(f"Running STAR alignment for paired sample: {sample_id}")
            try:
                cmd = [
                    "conda", "run", "-n", conda_env,
                    "STAR", "--genomeDir", index_dir,
                    "--readFilesIn", r1_file, r2_file,
                    "--outFileNamePrefix", f"{sample_output_dir}/",
                    "--runThreadN", str(threads),
                    "--outSAMtype", "BAM", "SortedByCoordinate",
                    "--outSAMattributes", "Standard"
                ]
                
                # Handle gzipped files
                if r1_file.endswith(".gz"):
                    cmd.extend(["--readFilesCommand", "zcat"])
                
                subprocess.run(cmd, check=True)
                
                # Index the BAM file
                bam_file = os.path.join(sample_output_dir, "Aligned.sortedByCoord.out.bam")
                subprocess.run([
                    "conda", "run", "-n", conda_env, 
                    "samtools", "index", bam_file
                ], check=True)
                
                logger.info(f"Alignment completed for {sample_id}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to run alignment for {sample_id}: {e}")
                return False
    
    return True

def main():
    """Main function to run the STAR alignment workflow."""
    args = parse_arguments()
    
    # Validate inputs
    if not validate_inputs(args.genome, args.fastq_dir):
        sys.exit(1)
    
    # Setup conda environment
    if not setup_conda_environment(args.conda_env):
        sys.exit(1)
    
    # Get strandedness
    strandedness = get_strandedness()
    logger.info(f"Processing {strandedness}-end data")
    
    # Set up STAR index directory
    if args.star_index_dir:
        index_dir = args.star_index_dir
    else:
        genome_dir = os.path.dirname(os.path.abspath(args.genome))
        index_dir = os.path.join(genome_dir, "star_index")
    
    # Check if index exists and build if necessary
    index_exists = os.path.exists(index_dir) and os.path.exists(os.path.join(index_dir, "Genome"))
    
    if not index_exists or args.force_index:
        if not build_star_index(args.genome, index_dir, args.threads, args.conda_env):
            sys.exit(1)
    else:
        logger.info(f"Using existing STAR index at {index_dir}")
    
    # Get fastq files
    fastq_files = find_fastq_files(args.fastq_dir)
    logger.info(f"Found {len(fastq_files)} FASTQ files")
    
    # Run alignment
    if not run_star_alignment(index_dir, fastq_files, args.output_dir, 
                             strandedness, args.threads, args.conda_env):
        sys.exit(1)
    
    logger.info("STAR alignment completed successfully!")

if __name__ == "__main__":
    main()