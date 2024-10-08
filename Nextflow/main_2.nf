#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define your parameters and file paths
params.final_genomes = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/final_genomes"
params.pairs = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/batch1_pairs_old.csv"
params.bin_dir = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/bin"
params.alignment_output = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/alignments"
params.divergence_output = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/divergence"

process align_genomes {
    tag "$pair"
    clusterOptions '--time=1-00:00:00 -A highmem --ntasks 64 --mem-per-cpu=2G'

    input:
    tuple val(pair), path(genome1), path(genome2)

    output:
    tuple val(pair), path("${pair}.paf"), path("${pair}.srt.paf"), path("${pair}.var.txt")

    publishDir params.alignment_output, mode: 'copy'

    errorStrategy 'ignore' // Ignore errors and proceed with the workflow

    script:
    """
    # Run the alignment command and capture the exit status
    bash ${params.bin_dir}/align.sh ${pair} ${genome1} ${genome2}
    status=\$?

    # Check if the alignment was successful
    if [[ \$status -ne 0 ]]; then
        echo "Alignment failed for pair: ${pair} with exit code \$status" >> ${params.alignment_output}/failed_alignments.log

        # Create empty or placeholder files to avoid pipeline failure
        [[ ! -f ${pair}.paf ]] && touch ${pair}.paf
        [[ ! -f ${pair}.srt.paf ]] && touch ${pair}.srt.paf
        [[ ! -f ${pair}.var.txt ]] && touch ${pair}.var.txt
    fi
    """
}

process calculate_divergence {
    tag "$pair"
    clusterOptions '--ntasks 16 --time 1-00:00:00 -A johnwayne'

    input:
    tuple val(pair), path(srt_paf), path(var_file)

    output:
    path("${pair}_distance.csv")

    publishDir params.divergence_output, mode: 'copy'

    script:
    """
    echo "Processing pair: ${pair}"
    echo "PAF file: ${srt_paf}"
    echo "VAR file: ${var_file}"
    
    python3 ${params.bin_dir}/calc_div.py ${var_file} ${srt_paf} \$PWD ${pair}
    
    if [[ -f \$PWD/${pair}_distance.csv ]]; then
        echo "Output file ${pair}_distance.csv generated successfully."
    else
        echo "Error: Output file ${pair}_distance.csv not found."
        exit 1
    fi
    """
}

// Workflow for alignment and divergence calculation
workflow {
    // Process genome pairs from the pairs.csv file
    Channel.fromPath(params.pairs)
        | splitCsv(header: true)
        | map { row -> 
            def genome1 = file("${params.final_genomes}/${row.genome1}_final.fna")
            def genome2 = file("${params.final_genomes}/${row.genome2}_final.fna")
            [row.pair, genome1, genome2]
        }
        | filter { pair, genome1, genome2 -> genome1.exists() && genome2.exists() }
        | align_genomes
        | calculate_divergence
        | collect
}
