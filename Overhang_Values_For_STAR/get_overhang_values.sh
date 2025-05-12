#!/bin/bash

echo "-------------------------------------------------"
echo "07/28/2023"
echo "get_overhang_values.sh"
echo "Made by Gabriel Dall'Alba"
echo "Script to get overhang values required"
echo "for overhang parameter in STAR genome indexing"
echo "please provide folder containing fastq.gz files"
echo "-------------------------------------------------"
# directory containing fastq.gz files
dir=$1

# arrays to store all lengths and individual file averages
all_lengths=()
averages=()

# iterate over each fastq.gz file
for file in $dir/*.fastq.gz
do
  echo "Processing $file..."
  
  # use awk to process the sequence lines (every 4th line starting from the 2nd) and get lengths
  lengths=$(zcat $file | awk '{if(NR%4==2) print length($1)}')
  
  # calculate and print average for this file
  file_avg=$(echo $lengths | awk '{ sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF }')
  echo "Average read length for $file: $file_avg"
  averages+=($file_avg)
  
  # add lengths to all lengths array
  all_lengths+=($lengths)
  
  # print read length distribution for this file
  echo "Read length distribution for $file:"
  printf '%s\n' "${lengths[@]}" | sort -n | uniq -c | sort -bgr
done

# calculate and print overall average
overall_avg=$(echo ${all_lengths[*]} | awk '{ sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF }')
echo "Average read length over all files: $overall_avg"

# print read length distribution over all files
echo "Read length distribution over all files:"
printf '%s\n' "${all_lengths[@]}" | sort -n | uniq -c | sort -bgr

# print averages for each file
echo "Averages for each file:"
printf '%s\n' "${averages[@]}"
