#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.csvA = "/scratch/negishi/allen715/nextflow_templates/gd_pipe/genomes.csv"
params.temp_dir = "/scratch/negishi/allen715/nextflow_templates/gd_pipe/tmp"
params.genomes = "/scratch/negishi/allen715/nextflow_templates/gd_pipe/genomes"
params.final_genomes = "/scratch/negishi/allen715/nextflow_templates/gd_pipe/final_genomes"
params.pairs = "/scratch/negishi/allen715/nextflow_templates/gd_pipe/pairs.csv"
params.bin_dir = "/scratch/negishi/allen715/nextflow_templates/gd_pipe/bin"
params.alignment_output = "/scratch/negishi/allen715/nextflow_templates/gd_pipe/alignments"
params.divergence_output = "/scratch/negishi/allen715/nextflow_templates/gd_pipe/divergence"

process get_genomes {
    tag "$accession"
    clusterOptions '--ntasks 16 --time 01:00:00 -A johnwayne'

    input:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask)

    output:
    tuple val(accession), val(genus), val(species), val(assembly), val(full_scientific_name), val(taxonomic_group), val(mask), path("${genus}_${species}_${accession}.fna")

    script:
    """
    sh ${params.bin_dir}/get_genomes_nf.sh $params.temp_dir \$PWD $accession $genus $species $assembly
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

    publishDir params.final_genomes, mode: 'copy'

    script:
    """
    cp ${masked_file} ${genus}_${species}_${accession}.final.fna
    mv ${genus}_${species}_${accession}.final.fna ${accession}_final.fna
    """
}

process align_genomes {
    tag "$pair"
    clusterOptions '--ntasks 64 --time 1-00:00:00 -A highmem'

    input:
    tuple val(pair), path(genome1), path(genome2)

    output:
    tuple val(pair), path("${pair}.paf"), path("${pair}.srt.paf"), path("${pair}.var.txt")

    publishDir params.alignment_output, mode: 'copy'

    script:
    """
    bash ${params.bin_dir}/align.sh ${pair} ${genome1} ${genome2}
    """
}

process calculate_divergence {
    tag "$pair"
    clusterOptions '--ntasks 16 --time 1-00:00:00 -A highmem'

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
    
    # Check if the output file is generated
    if [[ -f \$PWD/${pair}_distance.csv ]]; then
        echo "Output file ${pair}_distance.csv generated successfully."
    else
        echo "Error: Output file ${pair}_distance.csv not found."
        exit 1
    fi
    
    # Print the content of the generated file for debugging
    echo "Content of ${pair}_distance.csv:"
    cat \$PWD/${pair}_distance.csv
    
    """
}

process concatenate_divergence {
    clusterOptions '--ntasks 16 --time 01:00:00 -A standby'

    input:
    path divergence_files

    output:
    path("${params.divergence_output}/all_divergence.csv")

    script:
    """
    echo "Files to be processed:"
    ls ${divergence_files}

    # Ensure the output directory exists
    mkdir -p ${params.divergence_output}

    echo "pair,raw,k2p,k3p,jc" > ${params.divergence_output}/all_divergence.csv
    for file in ${divergence_files}; do
        echo "Processing file: \$file"
        pair=\$(basename \${file} _distance.csv)
        distances=\$(tail -n +2 \${file} | tr '\\n' ',' | sed 's/,\$//')
        echo "\${pair},\${distances}" >> ${params.divergence_output}/all_divergence.csv
    done
    """
}

workflow {
    print "Welcome to the workflow operating on an input CSV"
    print 'You will be working on the following file: ' + params.csvA

    Channel.fromPath(params.csvA) \
        | splitCsv(header: true) \
        | map { row -> [row.accession, row.genus, row.species, row.assembly, row.full_scientific_name, row.taxonomic_group, row.mask] } \
        | set { csvA }

    get_genomes(csvA) \
        | mask \
        | move_to_final \
        | set { final_genomes }

    Channel.fromPath(params.pairs) \
    | splitCsv(header: true) \
    | map { row -> 
        def genome1 = file("${params.final_genomes}/${row.genome1}")
        def genome2 = file("${params.final_genomes}/${row.genome2}")
        println "Processing pair: ${row.pair}, Genome1: ${genome1}, Genome2: ${genome2}"
        [row.pair, genome1, genome2]
    } \
    | view { pair, genome1, genome2 -> 
        println "Aligning genomes for pair: ${pair}"
        println "Full path for Genome 1: ${genome1.toAbsolutePath()}"
        println "Full path for Genome 2: ${genome2.toAbsolutePath()}"
        return [pair, genome1, genome2]
    } \
    | align_genomes \
    | set { aligned_pairs }

    aligned_pairs \
        | map { pair, paf, srt_paf, var_file -> 
            [pair, srt_paf, var_file.toString().replaceFirst('.srt.paf', '.var.txt')]
        } \
        | calculate_divergence \
        | collect() \
        | view { divergence_files -> 
            println "Divergence files passed to concatenate_divergence: ${divergence_files}"
        } \
        | concatenate_divergence
}

