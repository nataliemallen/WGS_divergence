#!/usr/bin/python3
import csv
import math
import os
import sys

def calculate_distances(var_file, paf_file, output_dir, pair):
    transition_count = 0
    transversion_count = 0
    total_differences = 0
    total_sites = 0

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

    p = transition_count / total_sites
    q = transversion_count / total_sites
    r = total_differences / total_sites

    raw_distance = r
    k2p_distance = (-0.5 * math.log(1 - (2 * p - q))) - (0.25 * math.log(1 - (2 * q)))
    k3p_distance = (-0.5 * math.log(1 - (2 * p - q))) - (0.20 * math.log(1 - (2 * q)))
    jc_distance = -0.75 * math.log(1 - (4/3) * r)

    output_file = os.path.join(output_dir, f"{pair}_distance.csv")
    with open(output_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["pair", "raw", "k2p", "k3p", "jc"])
        writer.writerow([pair, raw_distance, k2p_distance, k3p_distance, jc_distance])
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    var_file = sys.argv[1]
    paf_file = sys.argv[2]
    output_dir = sys.argv[3]
    pair = sys.argv[4]
    calculate_distances(var_file, paf_file, output_dir, pair)
