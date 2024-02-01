import csv
import math
import os
import sys

def calculate_distances(var_file, paf_file):
    # Initialize counts for transitions, transversions, and total differences
    transition_count = 0
    transversion_count = 0
    total_differences = 0
    total_sites = 0

    # Define sets for transition and transversion bases
    transitions = {"ag", "ga", "ct", "tc"}
    transversions = {"ac", "ca", "at", "ta", "gc", "cg", "gt", "tg"}

    with open(var_file, 'r') as file:
        for line in file:
            if not line.startswith('V'):
                continue

            columns = line.strip().split()

            if columns[4] != '1':
                continue

            ref_allele = columns[6].lower()
            alt_allele = columns[7].lower()

            if len(ref_allele) != 1 or len(alt_allele) != 1:
                continue

            change = ref_allele + alt_allele
            if change in transitions:
                transition_count += 1
            elif change in transversions:
                transversion_count += 1

            total_differences += 1

    with open(paf_file, 'r') as file:
        for line in file:
            columns = line.strip().split()
            start = int(columns[7])
            end = int(columns[8])
            total_sites += end - start

    # Calculate the proportions of transitions, transversions, and total differences
    p = transition_count / total_sites
    q = transversion_count / total_sites
    r = total_differences / total_sites

    # Calculate the distances
    raw_distance = r
    k2p_distance = (-0.5 * math.log(1 - (2 * p - q))) - (0.25 * math.log(1 - (2 * q)))
    k3p_distance = (-0.5 * math.log(1 - (2 * p - q))) - (0.20 * math.log(1 - (2 * q)))
    jc_distance = -0.75 * math.log(1 - (4/3) * r)

    # Write the distances to the output CSV file
    output_file = os.path.splitext(var_file)[0] + '_distance.csv'
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Raw Distance', raw_distance])
        writer.writerow(['K2P Distance', k2p_distance])
        writer.writerow(['K3P Distance', k3p_distance])
        writer.writerow(['JC Distance', jc_distance])

# Call the function with your input and output file paths
#calculate_distances('falcons.var.txt', 'falcons.srt.paf')

if __name__ == "__main__":
    # Check if the correct number of command line arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python3 script.py <input1> <input2>")
        sys.exit(1)

    # Get input file paths from command line arguments
    input_file1 = sys.argv[1]
    input_file2 = sys.argv[2]

    # Call the function with the input file paths
    calculate_distances(input_file1, input_file2)
