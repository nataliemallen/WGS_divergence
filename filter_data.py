#filters ncbi genome spreadsheet to contain only genera with >1 species

import csv
import os
from collections import defaultdict

os.getcwd()
os.chdir('/Users/natal/Documents/Purdue/Genomic divergence/dataset')

#input and output files
input_file = 'ncbi_dataset_vertebrates.tsv'
output_file = 'vertebrates_congeneric.tsv'

#dictionar for genus counts
genus_counts = defaultdict(int)

#read input and count genus occurences
with open(input_file, 'r', newline='') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    for row in reader:
        genus_counts[row['Genus']] += 1

#write rows with Genus occurring more than once to output file
with open(output_file, 'w', newline='') as f_out:
    writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter='\t')
    writer.writeheader()
    
    with open(input_file, 'r', newline='') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        for row in reader:
            if genus_counts[row['Genus']] > 1:
                writer.writerow(row)

print(f"Filtered data written to '{output_file}'.")
