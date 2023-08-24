# WGS_divergence
Scripts for calculating genomic divergence between pairs of congeneric species
General steps:
1. Obtain a csv file of all the genomes you want, including genus, species, genbank ID, and assembly number
2. Use NCBI dataseets to download all genomes and sort them into folders by genus
3. Use python to get a cvs with all possible species pairs within a genus 
4. For each pair, pick one genome and split into 1000 bp segments
5. 
