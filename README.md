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
  "conflicting" variant pairs produced with the
`get_conflicting_variants_indices.pl` script to filter phased haplotype blocks
assembled with HapCUT2. It also carries out filtering of haplotype blocks based
on HapCUT2 switch/mismatch qualities.



* the `find_haplotypic.counterparts.pl` Perl script takes as input reconstructed
 sequences of haplotypes of several individuals from the same locus and
identifies cases where both haplotypes of the same individual have
reciprocal closest counterparts in other individuals (without 
examining bootstrap support).

* the `print_nodes_2leaves.py` Python script reads in a user-specified
  phylogenetic tree in the Newick format and searches for nodes with two
leaves. For each such node, prints out its leaves and the corresponding
bootstrap support. Can be used in conjunction with the 
`find_haplotypic.counterparts.pl` script to find groupings of haplotypes
with decent bootstrap support.

* the `compute_H_scores.pl` Perl script searches an input VCF file for
 triallelic sites carrying all three heterozygous genotypes that harbor 
only one private heterozygous genotype and computes H-scores using this
subset of sites. 


# Comments

This repository contains stand-alone scripts and no installation is needed.
For all scripts, the expected run time for demo on a "normal" desktop computer
is less than 60 seconds.

Detailed information on the data requirements and dependencies for each
script is provided in the `README.md` located in the corresponding
subdirectory.

The provided scripts can be run on the actual data in the same way as they
are run on the sample input data (which represent subsets of real data
used in the study).
If you have any questions or require further information, please, contact 
the author (Olga Vakhrusheva, vakh57@gmail.com). 

