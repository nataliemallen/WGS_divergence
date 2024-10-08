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
