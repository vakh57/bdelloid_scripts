# Description:

The `compute_genotypic_distances.pl` script computes pairwise genotypic
distances for all pairs of samples present in the input VCF file.

If an input VCF file contains whole-genome calls, pairwise whole-genome
distances will be computed.  Alternatively, distances between samples at
individual loci can be computed: in this case an input VCF should contain only
entries corresponding to the locus of interest.

The distance between a pair of individuals is calculated as a sum of genotypic
distances for each polymorphic site present in the input VCF file normalized by
2\*n, where n corresponds to the total number of genotyped sites in the VCF file
(both polymorphic and monomorphic). The pairwise genotypic distance at a single
polymorphic site is defined as the difference between the harbored numbers of
non-reference variants. For example, the distance between the genotype 0/1 and
the genotype 0/0 is 1, and the distance between genotypes 0/0 and 1/1 is 2.


# Requirements:

* This script requires Perl 5 (tested with Perl v5.10.1 under Scientific Linux
  release 6.8, Carbon).
* This script is supposed to work only with biallelic SNP sites and monomorphic
  sites genotyped in all the samples contained in the VCF file. 
* Multiallelic sites and INDEL sites are not supposed to be present in the
  input VCF file.
* Only the following genotypes are allowed:  0/0 0/1 1/1.
* The input VCF file must not contain missing genotypes.
* Only a single entry is supposed to be present for each included genomic site.
* Sample names must not contain dots (.) 


# Usage:

`perl compute_genotypic_distances.pl input_vcf`

where: 
`input_vcf` is the name of the input VCF file


# Output files:

The script produces three output files:

* output file with the extension `.pair.gntp.frac.allsites.dist.txt` contains
  pairwise genotypic distances normalized by the total effective number of
analyzed sites (both polymorphic and monomorphic). If the input VCF file meets
all above-mentioned requirements, the total effective number of analyzed sites
is equal to 2\*n, where n is the number of sites present in the VCF file.

* output file with the extension `.pair.gntp.abs.dist.txt` contains absolute
  pairwise genotypic distances (not normalized).

* output file with the extension `.pair.gntp.frac.segsites.dist.txt` contains
  pairwise genotypic distances normalized by the effective number of analyzed
polymorphic sites. If the input VCF file meets all above-mentioned
requirements, the effective number of analyzed polymorphic sites is equal to
2\*n\_var, where n\_var is the number of biallelic polymorphic sites present in
the VCF file.


# Sample input data:

A sample input VCF file is included (`sample.vcf`) along with the script.

The command to run the script with the sample VCF file:
`perl compute_genotypic_distances.pl sample.vcf`

The corresponding sample output files are provided in the subdirectory
`sample_output`.
