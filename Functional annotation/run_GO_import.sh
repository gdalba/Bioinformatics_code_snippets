#!/bin/bash

# Directory containing the .xlsx files
DIRECTORY="$1"

# Activate Conda environment
echo "Activating Conda environment..."
source /project/(REMOVED_FOR_PRIVACY)/gdalba/anaconda3/etc/profile.d/conda.sh
conda activate functional

if [ $? -eq 0 ]; then
    echo "Environment activated successfully."
else
    echo "Failed to activate environment."
    exit 1
fi

for file in "$DIRECTORY"/*.xlsx
do
    echo "Processing file: $file"
    python import_GOs_to_table.py "$file"
    
    if [ $? -eq 0 ]; then
        echo "Script completed successfully for $file."
    else
        echo "Script encountered an error processing $file."
        continue
    fi
done

echo "All files processed successfully."
