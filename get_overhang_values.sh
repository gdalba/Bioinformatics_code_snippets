#!/bin/bash

echo "-------------------------------------------------"
echo "ctenotation v1 - 07/28/2023"
echo "get_overhang_values.sh"
echo "Made by Gabriel Dall'Alba"
echo "Script to get overhang values required"
echo "for overhang parameter in STAR genome indexing"
echo "please provide folder containing fastq.gz files"
echo "-------------------------------------------------"

# directory containing fastq.gz files
dir=$1

# output directory
out_dir=$2
mkdir -p $out_dir

# arrays to store all lengths and individual file averages
all_lengths=()
averages=()

# iterate over each fastq.gz file
for file in $dir/*.fastq.gz
do
  # get the filename without extension and path
  filename=$(basename -- "$file")
  filename="${filename%.*.*}"
  
  echo "Processing $file..."
  
  # use awk to process the sequence lines (every 4th line starting from the 2nd) and get lengths
  lengths=$(zcat $file | awk '{if(NR%4==2) print length($1)}')
  
  # calculate and print average for this file
  file_avg=$(echo $lengths | awk '{ sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF }')
  echo "Average read length for $file: $file_avg"
  averages+=($file_avg)
  
  # add lengths to all lengths array
  all_lengths+=($lengths)
  
  # save average and read length distribution for this file to CSV
  echo "Read length,Count" > "$out_dir/${filename}_distribution.csv"
  printf '%s\n' "${lengths[@]}" | sort -n | uniq -c | awk '{print $2","$1}' | sort -bgr >> "$out_dir/${filename}_distribution.csv"
  echo "Filename,Average" > "$out_dir/${filename}_average.csv"
  echo "${filename},${file_avg}" >> "$out_dir/${filename}_average.csv"
done

# calculate and print overall average
overall_avg=$(echo ${all_lengths[*]} | awk '{ sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF }')
echo "Average read length over all files: $overall_avg"

# save overall average and read length distribution to CSV
echo "Read length,Count" > "$out_dir/overall_distribution.csv"
printf '%s\n' "${all_lengths[@]}" | sort -n | uniq -c | awk '{print $2","$1}' | sort -bgr >> "$out_dir/overall_distribution.csv"
echo "Filename,Average" > "$out_dir/overall_average.csv"
echo "Overall,${overall_avg}" >> "$out_dir/overall_average.csv"

# print averages for each file
echo "Averages for each file:"
printf '%s\n' "${averages[@]}"