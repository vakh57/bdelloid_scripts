# Description:

The `compute_H_scores.pl` searches an input VCF file for triallelic sites 
carrying all three heterozygous genotypes that harbor only one private 
heterozygous genotype. That is, such sites where three heterozygotes are found
among the analyzed individuals and the least frequent of the three heterozygous 
genotypes is present in a single individual with the next frequent genotype 
present at least in two individuals. 
This subset of sites is used to compute H-scores. 
The logic and procedure behind computing H-scores is detailed in Supplementary 
Information. 

# Requirements:

* This script requires Perl 5 (tested with Perl v5.10.1 under Scientific Linux
  release 6.8, Carbon).
* [Data::Dumper](https://perldoc.perl.org/Data/Dumper.html) Perl module
* Input VCF file is supposed to contain only triallelic sites carrying three 
  different heterozygous genotypes. Other types of sites will be skipped, if present.
* Three different heterozygous genotypes have to be encoded as '0/1', '0/2' and '1/2'.
* This script is supposed to work only with sites genotyped in all the samples 
  contained in the VCF file. 
* The input VCF file must not contain missing genotypes.
* Only the following genotypes are allowed:  '0/0', '0/1', '1/1', '0/2', '2/2', '1/2'.


# Usage:

`perl compute_H_scores.pl input_3het_vcf`

where: 
`input_3het_vcf` is the name of the input VCF file containing triallelic sites
carrying all three possible heterozygous genotypes ('0/1', '0/2' and '1/2') 
among the analyzed individuals.


# Output files:

The script produces two output files:

* output file with the extension `.private.het.counts.txt` 
 contains per sample numbers of sites with three heterozygotes where 
 the least frequent heterozygous genotype is private to this sample.
 Only those sites with three heterozygotes where the least frequent 
 of the three heterozygous genotypes is present in a single individual 
 with the next frequent genotype present at least in two individuals
 are cosidered.

* output file with the extension `.H.scores.txt` contains 
 H-scores computed for all possible pairs of samples.


# Sample input data:

A sample input VCF file containing triallelic SNPs represented by
three heterozygous genotypes among the samples present in the
VCF file is included (`sample.3.hets.vcf`) along with the script.

The command to run the script with the sample VCF file:
`perl compute_H_scores.pl sample.3.hets.vcf`

The corresponding sample output files are provided in the subdirectory
`sample_output`.
