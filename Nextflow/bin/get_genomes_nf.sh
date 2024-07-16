#!/bin/bash
module load anaconda/2022.10-py39
#must install ncbi_datasets using conda, instructions here https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
conda activate ncbi_datasets
#on my home directory I did the following per the instructions at the link above
#
#First, create a conda environment: conda create -n ncbi_datasets
#Then, activate your new environment: conda activate ncbi_datasets
#Finally, install the datasets conda package: conda install -c conda-forge ncbi-datasets-cli
#nextflow works after doing the installation manually once in home directory

TEMP_DIR=$1
WORK_DIR=$2
ACCESSION=$3 
GENUS=$4
SPECIES=$5
ASSEMBLY=$6

cd $TEMP_DIR/
mkdir -p ${ACCESSION}
cd $ACCESSION/

# Download the genome zip file and unzip it
datasets download genome accession "$ACCESSION" --filename "${GENUS}_${SPECIES}_${ACCESSION}.zip"
unzip ${GENUS}_${SPECIES}_${ACCESSION}.zip

# Copy the genome file to the work directory and rename it
cp "$TEMP_DIR/${ACCESSION}/ncbi_dataset/data/$ACCESSION/${ACCESSION}_${ASSEMBLY}_genomic.fna" "$WORK_DIR/${GENUS}_${SPECIES}_${ACCESSION}.fna"

# Clear unneeded files out of temp directory
rm -r "$TEMP_DIR/${ACCESSION}"
