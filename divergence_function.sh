#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=divergence
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module load biocontainers
module load biopython
module load r/4.3.1

run_processing() {
  NAME="$1"

  for NUM in 5 10 20; do
    NAME_var_NUM="/scratch/negishi/allen715/GD_tests/minimap/tests/${NAME}_${NUM}.var.txt"
    NAME_paf_NUM="/scratch/negishi/allen715/GD_tests/minimap/tests/${NAME}_${NUM}.srt.paf"

    python calc_all.py "$NAME_var_NUM" "$NAME_paf_NUM"
    Rscript paf_colinear_dotplot2.R -i "$NAME_paf_NUM" -o out -s -t -m 500 -q 300000 -l
  done
}

# Example usage:
#run_processing "chicken_grouse_masked"
