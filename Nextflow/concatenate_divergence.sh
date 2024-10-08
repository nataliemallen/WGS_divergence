#!/bin/bash

# Define the directory containing the distance files
distance_dir="/path/to/your/distance/files"
output_file="all_distances.csv"

# Write the header to the output file
echo "pair,raw,k2p,k3p,jc" > $output_file

# Append the contents of each distance.csv file to the output file (skipping the header)
for file in "$distance_dir"/*_distance.csv; do
    tail -n +2 "$file" >> $output_file
done

echo "Concatenation complete. Output written to $output_file"
