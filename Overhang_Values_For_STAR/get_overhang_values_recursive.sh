#!/bin/bash

# Script path
script_path=$1

# Shift the first argument
shift 1

# Iterate over directories
for dir in "$@"; do
  echo "Processing $dir..."
  
  # Create output directory named after the input directory
  out_dir="output_$(basename $dir)"
  mkdir -p "$out_dir"
  
  $script_path "$dir" "$out_dir"
done
