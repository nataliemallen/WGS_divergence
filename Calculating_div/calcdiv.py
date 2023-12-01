import csv
import math

def calculate_k2p_distance(var_file, output_file):
    # Initialize counts for transitions and transversions
    transition_count = 0
    transversion_count = 0

    # Define sets for transition and transversion bases
    transitions = {"ag", "ga", "ct", "tc"}
    transversions = {"ac", "ca", "at", "ta", "gc", "cg", "gt", "tg"}

    with open(var_file, 'r') as file:
        for line in file:
            # Skip lines that do not start with 'V'
            if not line.startswith('V'):
                continue

            # Split the line into columns
            columns = line.strip().split()

            # Skip lines where column 5 is not one
            if columns[4] != '1':
                continue

            # Get the REF and ALT alleles
            ref_allele = columns[6].lower()
            alt_allele = columns[7].lower()

            # Skip lines where REF or ALT is not a single base
            if len(ref_allele) != 1 or len(alt_allele) != 1:
                continue

            # Check if the change is a transition or transversion
            change = ref_allele + alt_allele
            if change in transitions:
                transition_count += 1
            elif change in transversions:
                transversion_count += 1

    # Calculate the proportions of transitions and transversions
    p = transition_count / (transition_count + transversion_count)
    q = transversion_count / (transition_count + transversion_count)

    # Print the counts and proportions for debugging
    print(f'Transitions: {transition_count}, Transversions: {transversion_count}, p: {p}, q: {q}')

    # Check if the argument of the math.log function is greater than zero
    if (1 - 2*p - q) * math.sqrt(1 - 2*q) > 0:
        # Calculate the K2P distance
        k2p_distance = -0.5 * math.log((1 - 2*p - q) * math.sqrt(1 - 2*q))
    else:
        k2p_distance = 'undefined'

    # Write the K2P distance to the output CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['K2P Distance'])
        writer.writerow([k2p_distance])

# Call the function with your input and output file paths
calculate_k2p_distance('falcons.var.txt', 'output.csv')
