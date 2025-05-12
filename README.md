# Bioinformatics Code Snippets

This repository contains code that was developed by me over the last few years. Not everything is perfect nor For the Bioinformatics Scientist job.

---
## Data Processing

On-demand small code to address specific needs. Usually not elegant.

`DOWNLOAD_ACCESSIONS.sh`: For download of multiple NCBI datasets in one go.

`bol_gff_cleaner_v1.py`: Formatting of .gff3 file with a missing column.

`ctenotation_gff_statistics.py`: Compiles statistics of .gff3 file. Creates a database .db file with gffutils for faster access. 

`extract_and_get_fasta.py`: If a gene had no annotation, I'd get the sequences to try different approaches later on.

---
## Functional Annotation

After compiling functional annotation about gene models, these were useful to translate IDs to full-length names. (They suffer from lazyness of removing unnecesssary parts)

`GOs_by_obo_optimized.py`: Makes use of GO's .obo file to access entries locally instead of requesting through slow API.

`KEGG_IDs.py` & `KEGG_PATHWAYS.py`: Same goal as above but through API (using bioservices).

To process hundreds of thousands of entries through API, I deploy `import_GOs_to_Table.py` called through `run_GO_import.sh`

---
## Overhang_Values_For_STAR
	
`get_overhang_values.sh` and `get_overhang_values_save_on_tables.sh` major difference is the latter will save outputs to spreadsheets. It also does some things better than the other script, such as handling of files. I needed this to manually modify parameters within STAR aligner. These tools process multiple FASTQ files (compressed .fastq.gz), calculate average read length for each file, and provides read length distributions to identify dominant lengths.

Usage:

```bash
./get_overhang_values_save_on_tables.sh /path/to/fastq/directory /path/to/output/directory
```

Finally, this little gem (`get_overhang_values_recursive.sh`) is a clever wrapper that adds batch processing capabilities to the original read length analysis tools. This script takes the original overhang calculation script and applies it to multiple directories in one go.

Usage:

```bash
./get_overhang_values_recursive.sh ./get_overhang_values_save_on_tables.sh /path/to/dataset1 /path/to/dataset2 /path/to/dataset3
```

---

## PBS

Contain my modularization efforts to build `.pbs` files and automatically submit them to UBC's SOCKEYE HPC cluster. `PROJECT_ENV(REMOVED_FOR_PRIVACY)` is intentionally modified to not give sensitive information about my lab's project name.

---

## Visualization_tools

Case-specific needs of generating plots.

`ORIGINAL_TE_log_scatterplot.py`: For pairwise comparisons across genomes. At the time I wanted to compare three genomes pairwise. Output is under embargo and cannot be displayed.

---

## Misc

Contains non-categorize examples of code.

---
 
