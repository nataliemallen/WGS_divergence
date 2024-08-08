#!/bin/bash

module load biocontainers
module load minimap2

# Function to run minimap2 and handle errors
run_minimap() {
    if [ "$#" -ne 3 ]; then
        echo "Usage: run_minimap <pair> <file1> <file2>"
        return 1
    fi

    pair="$1"
    file1="$2"
    file2="$3"

    output_dir="/scratch/negishi/allen715/nextflow_templates/gd_pipe/alignments"
    mkdir -p "$output_dir"

    echo "Running minimap2 for $file1 and $file2"
    minimap2 -cx asm10 -t64 --cs "$file1" "$file2" > "${output_dir}/${pair}.paf"
    if [ $? -ne 0 ]; then
        echo "Error: minimap2 failed"
        return 1
    fi

    echo "Sorting PAF file"
    sort -k6,6 -k8,8n "${output_dir}/${pair}.paf" > "${output_dir}/${pair}.srt.paf"
    if [ $? -ne 0 ]; then
        echo "Error: sorting PAF file failed"
        return 1
    fi

    echo "Calling variants with paftools.js"
    paftools.js call "${output_dir}/${pair}.srt.paf" > "${output_dir}/${pair}.var.txt"
    if [ $? -ne 0 ]; then
        echo "Error: paftools.js call failed"
        return 1
    fi

    echo "Minimap2 alignment and variant calling completed successfully"

    # Move output files to current working directory
    cp "${output_dir}/${pair}.paf" .
    cp "${output_dir}/${pair}.srt.paf" .
    cp "${output_dir}/${pair}.var.txt" .
}

# Call the function with the provided arguments
run_minimap "$1" "$2" "$3"
