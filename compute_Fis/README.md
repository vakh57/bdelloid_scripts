# Description:

The `compute_Fis.sh` script computes values of Fis statistics for individual 
biallelic sites based on observed and expected numbers of heterozygous genotypes. 
The script uses data on observed and expected proportions of heterozygotes 
among the individuals of interest (typically, drawn from the same population). 
An input file must be in the format produced by the `populations` 
[program](http://catchenlab.life.illinois.edu/stacks/comp/populations.php), 
a part of the `Stacks` software pipeline. 

For each biallelic site present in the input file, the script computes the 
corresponding value of the inbreeding coefficient Fis.

Fis is computed as 1-Ho/He, where Ho and He stand for the observed and expected 
numbers of heterozygous genotypes respectively. 
Fractional expected genotype counts are rounded to the nearest integer number. 

This procedure was designed for small sample sizes, where fractional expected 
genotype counts can significantly bias Fis estimates. 
For example, with the sample size of 8 individuals and a minor allele count 
of 6 (out of 16 alleles), the unrounded value of He is 3.75, which of course 
can never be observed. 
To avoid artifacts in Fis calculation that can stem from using raw He values, 
the script uses rounded values of He.


# Requirements:

* This script requires bash, cut, sed, grep and Awk (tested with bash 4.1.2, 
  cut 8.4, sed 4.2.1, grep 2.20 and Awk 3.1.7 under Scientific Linux release 6.8, 
  Carbon).
* Input file containing data on observed and expected proportions of heterozygous
  genotypes for biallelic loci of interest created with the `populations` 
  [program](http://catchenlab.life.illinois.edu/stacks/comp/populations.php), 
  a part of the `Stacks` software pipeline.

  An input file for the script can be produced from a multi-sample VCF file 
  with the following command:
  `populations -V input_vcf -O output_dir --fstats`

  where:	       
  `input_vcf` and `output_dir` are the names of the input VCF file and 
  output folder respectively. Among files created by the populations program, 
  select the one with the extension `.p.sumstats.tsv`.
  The script has been tested on output files created with version 2.4 of 
  the `populations` program.


# Usage:

`bash compute_Fis.sh sample_tsv`

where:
`sample_tsv` is the name of the input file produced
with the `populations` program.

# Output files:

The script produces a single tab-delimited output file with the extension 
`.Fis.txt` containing computed Fis values.
Each line of the output file corresponds to a single biallelic locus 
listed in the input file.

## Format of output files:

The output file contains 18 columns. Fis values can be found in the last column. 
Remaining columns contain other related data, mostly from the input file created 
by `populations` and are specified in the output file header.

# Sample input data:

A sample input file (`sample.p.sumstats.tsv`) 
created by the `populations` program
is included along with the script.

The command to run the script with the sample input file:
`bash compute_Fis.sh sample.p.sumstats.tsv`

The corresponding sample output file is provided in the subdirectory
`sample_output`.

