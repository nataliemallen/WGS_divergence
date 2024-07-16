#!/bin/bash

module load biocontainers
module load biopython
module load repeatmasker

WORK_DIR=$1
GENOMES_DIR=$2
ACCESSION=$3
GENUS=$4
SPECIES=$5
ASSEMBLY=$6
TAX_GROUP=$7

# Determine the species for RepeatMasker
case "$TAX_GROUP" in
  bird)
    MASK_SPECIES="avian"
    ;;
  mammal)
    MASK_SPECIES="mammals"
    ;;
  fish|reptile|amphibian)
    MASK_SPECIES="vertebrates"
    ;;
  *)
    echo "Unknown taxonomic group: $TAX_GROUP"
    exit 1
    ;;
esac

# Define input and output files
INPUT_GENOME="$GENOMES_DIR/${GENUS}_${SPECIES}_${ACCESSION}.nm.fna"
OUTPUT_GENOME="${WORK_DIR}/${GENUS}_${SPECIES}_${ACCESSION}.masked.fna"

# Run RepeatMasker
RepeatMasker -pa 64 -qq -species "$MASK_SPECIES" -dir "$WORK_DIR" "$INPUT_GENOME"

# Move the masked genome to the final location
mv "${INPUT_GENOME}.masked" "$OUTPUT_GENOME"
mv "$OUTPUT_GENOME" "${GENOMES_DIR}/${GENUS}_${SPECIES}_${ACCESSION}.masked.fna"
