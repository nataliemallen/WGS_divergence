#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//paths to files
//for genomes.csv and pairs.csv, find and replace spaces with _

// Define your parameters and file paths
params.csvA = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/batch1_genomes_old.csv"
params.temp_dir = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/tmp"
params.genomes = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/genomes"
params.final_genomes = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/final_genomes"
params.pairs = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/batch1_pairs_old.csv"
params.bin_dir = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/bin"
params.alignment_output = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/alignments"
params.divergence_output = "/scratch/negishi/allen715/shared/nextflow/gd_pipe/divergence"
params.error_log = "${params.final_genomes}/download_errors.log"

process get_genomes {
    tag "$accession"
    clusterOptions '--ntasks 16 --time 01:00:00 -A johnwayne'
    errorStrategy 'ignore'

    input:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask)

    output:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask), path("${genus}_${species}_${accession}.fna"), optional: true

    script:
    """
    bash ${params.bin_dir}/get_genomes_wrapper.sh "$params.temp_dir" "\$PWD" "$accession" "$genus" "$species" "$assembly" "$params.bin_dir/get_genomes_nf.sh" "$params.error_log"
    """
}

process remove_mito {
    tag "$accession"
    clusterOptions '--ntasks 16 --time 01:00:00 -A johnwayne'

    input:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask), path(genome_file)

    output:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask), path("${genome_file.baseName}.nm.fna")

    script:
    """
    python3 ${params.bin_dir}/remove_mito.py ${genome_file}
    """
}

process mask {
    tag "$accession"
    clusterOptions '--ntasks 16 --time 01:00:00 -A johnwayne'

    input:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask), path(filtered_genome)

    output:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask), path("${accession}_masked.fna")

    script:
    """
    set -x
    echo "Masking process started for: ${accession}"
    echo "genus: ${genus}, species: ${species}, assembly: ${assembly}, taxonomic_group: ${taxonomic_group}, mask: ${mask}"

    if [ "${mask}" == "1" ]
    then
        sh ${params.bin_dir}/mask.sh \$PWD \$PWD ${accession} ${genus} ${species} ${assembly} ${taxonomic_group}
        mv ${genus}_${species}_${accession}.masked.fna ${accession}_masked.fna
        echo "File for ${accession} was masked" >> ${params.final_genomes}/log.txt
    else
        cp ${filtered_genome} ${accession}_masked.fna
        echo "File for ${accession} skipped the masking step" >> ${params.final_genomes}/log.txt
    fi
    """
}

process move_to_final {
    tag "$accession"
    clusterOptions '--ntasks 16 --time 01:00:00 -A johnwayne'

    input:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask), path(masked_file)

    output:
    tuple val(accession), path("${accession}_final.fna")

    publishDir params.final_genomes, mode: 'copy'  // or 'symlink' for linking instead of copying

    script:
    """
    # Instead of moving the file, simply rename it and let publishDir handle the copy
    cp "${masked_file}" "${accession}_final.fna"

    # Verify if the final file is created
    if [[ ! -f "${accession}_final.fna" ]]; then
        echo "Error: Final file ${accession}_final.fna was not created."
        exit 1
    fi
    """
}

workflow {
    // Fetch genomes and process them
    Channel.fromPath(params.csvA)
        | splitCsv(header: true)
        | map { row -> [row.accession, row.genus, row.species, row.assembly, row.full_scientific_name, row.taxonomic_group, row.mask] }
        | set { csvA }

    get_genomes(csvA)
        | filter { it.size() == 8 } // Ensure successful downloads
        | remove_mito
        | mask
        | move_to_final
        | set { final_genomes }
