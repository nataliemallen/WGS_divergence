If having problems with main.nf

Run main_1.nf first, then main_2.nf, then concatenate_divergence.sh

This will sometimes result in divergence values of 0 for any pairs where genomes are missing, alignment failed, or divergence could not be calculated. You may need to go back and check why divergence was not calculated for any such pairs
