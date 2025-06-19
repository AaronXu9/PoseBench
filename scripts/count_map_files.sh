#!/bin/bash

# Directory to check
base_dir="/home/aoxu/projects/PoseBench/forks/ICM/ICM_manual_docking_maps/plinder_set"

# Check if base directory exists
if [ ! -d "$base_dir" ]; then
    echo "Error: $base_dir directory not found"
    exit 1
fi

# Initialize counter
count=0

# Loop through all subdirectories
for subdir in "$base_dir"/*/; do
    # Check if any .map files exist in this subdirectory
    if ls "$subdir"*.map &> /dev/null; then
        count=$((count + 1))
        echo "Found .map files in: $subdir"
    fi
done

# Print the final count
echo "Total subdirectories containing .map files: $count"