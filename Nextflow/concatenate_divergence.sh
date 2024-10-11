# Write the header to the output file
echo "pair,raw,k2p,k3p,jc" > batch1_divergence_update.csv

# Append the contents of each distance.csv file to the output file (skipping the header)
for file in *_distance.csv; do
    tail -n +2 "$file" >>batch1_divergence_update.csv
done

echo "Concatenation complete. Output written to $output_file"
