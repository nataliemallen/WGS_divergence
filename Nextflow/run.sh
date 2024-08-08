#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH -t 14-00:00:00
#SBATCH --job-name=nextflow_job
#SBATCH --output=nextflow_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module load nextflow

nextflow run main.nf -c nextflow.config -resume

