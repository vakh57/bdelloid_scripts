# Introduction

This directory contains scripts that were used to process sequencing data for
*A. vaga* individuals L1-L11 analyzed in the
[manuscript](https://www.biorxiv.org/content/10.1101/489393v1). Each script is
in a separate subdirectory containing the script along with sample input and
output files.

# Available scripts

A list of provided scripts and short descriptions:

* the `compute_genotypic_distances.pl` Perl script computes pairwise genotypic
  distances for all pairs of samples present in a user-specified VCF file.

* the `compute_Fis.sh` Bash script calculates values of Fis statistics
  for biallelic loci. 

Two scripts to process and filter data on haplotype phasing generated with
[HapCUT2](https://github.com/vibansal/HapCUT2):
  
* the `get_conflicting_variants_indices.pl` Perl script looks for pairs of
  "conflicting" variants represented by more than two "haplotypes" in the
aligned reads from a single individual. This script reads in a user-specified
fragment file produced with the extractHAIRS command from HapCUT2.

* the `filter_hapcut2_haplotype_blocks.pl` Perl script uses a list of
  conflicting variant pairs produced with the
`get_conflicting_variants_indices.pl` script to filter phased haplotype blocks
assembled with HapCUT2. It also carries out filtering of haplotype blocks based
on HapCUT2 switch/mismatch qualities.


* the `print_nodes_2leaves.py` Python script reads in a user-specified
  phylogenetic tree in the Newick format and searches for nodes with two
leaves. For each such node, prints out its leaves and the corresponding
bootstrap support.

* the `compute_H_scores.pl` Perl script searches an input VCF file for
 triallelic sites carrying all three heterozygous genotypes that harbor 
only one private heterozygous genotype and computes H-scores using this
subset of sites. 


# Comments

Detailed information on the data requirements and dependencies for each
script is provided in the `README.md` located in the corresponding
subdirectory.