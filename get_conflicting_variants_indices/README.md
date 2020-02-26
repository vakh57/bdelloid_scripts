# Description:

The `get_conflicting_variants_indices.pl` script reads in a fragment file
produced with the extractHAIRS command from 
[HapCUT2](https://github.com/vibansal/HapCUT2).  For each pair of variants
simultaneously covered by at least one fragment, the script extracts
information on haplotypes supported by individual fragments.
The script looks for pairs of conflicting variants - that is, variants
represented by more than two haplotypes across different fragments.

A fragment corresponds to a single read or to a pair of reads generated from
both ends of the same DNA fragment.

The script reports two lists of indices corresponding to conflicting variant
pairs: 1) list of indices for all pairs of conflicting variants (irrespective
of read support)
2) list of indices for those pairs of conflicting variants that a) are present
in reads as three distinct haplotypes each supported by at least two fragments
or b) are present in reads as four haplotypes irrespective of the number of
fragments supporting different haplotypes.

**IMPORTANT:**
Index of a variant corresponds to its index number in the VCF file supplied to
extractHAIRS to produce the fragment file.
Index of a variant is NOT its chromosomal position!!!

For example, a SNP corresponding to the first record in a VCF file (i.e. first
non-header line) will get index 1, a SNP corresponding to the second record
will get index 2, etc.

# Requirements:

* This script requires Perl 5 (tested with Perl v5.10.1 under Scientific Linux
  release 6.8, Carbon).
* Input fragment matrix file generated with the extractHAIRS command from
  [HapCUT2](https://github.com/vibansal/HapCUT2).
* Input fragment matrix file is required to be generated with extractHAIRS
  based on a VCF file containing only biallelic single-nucleotide variants.    

# Usage:

`perl get_conflicting_variants_indices.pl input_fragment_matrix_file`

where:
`input_fragment_matrix_file` is the name of the input fragment file produced
with the `extractHAIRS` command from
[HapCUT2](https://github.com/vibansal/HapCUT2).

# Output files:

The script produces two output files:

* output file with the extension `.conflicting.indices.txt` contains a list of
  indices for all pairs of conflicting variants (irrespective of read support)
* output file with the extension `.conflicting.indices.supported.txt` contains
  a list of indices for those pairs of conflicting variants that a) are present
in reads as three distinct haplotypes each supported by at least two fragments
or b) are present in reads as four haplotypes irrespective of the number of
fragments supporting different haplotypes

## Format of output files:

Output files are semicolon separated.
Each line of the output file contains information on one pair of conflicting
variants.  Field 1 of each line is the index of the pair of conflicting
variants, constructed by concatenating the two corresponding variant indices
with a dot ('.') as a separator.  For example, record 186.194 corresponds to a
conflicting pair of variants with indices 186 and 194.

**IMPORTANT:**
Index of a variant corresponds to its index number in the VCF file supplied to
extractHAIRS to produce the fragment file.  Index of a variant is NOT its
chromosomal position!!!

Field 2 is the total number of haplotypes found in the reads for this pair of
conflicting variants.

Field 3 contains the list of haplotypes found in individual fragments for for
this pair of conflicting variants.  Haplotype encoding follows the VCF and
extractHAIRS encoding (0 stands for the reference and 1 for the alternative
allele). For example, haplotype encoding 01 corresponds to the case when the
haplotype carries the reference allele at the first variant in the pair and
alternative allele at the second variant.  For each haplotype listed in the
field 3, field 4 contains a number of fragments carrying this haplotype (the
order of haplotypes is the same as in field 3).

Field 5 gives the number of fragments supporting the least frequent haplotype
among those listed in field 3.

# Sample input data:

A sample input fragment matrix file is included
(`sample.hapcut2.fragment_matrix_file`) along with the script.

The command to run the script with the sample input fragment matrix file:
`perl get_conflicting_variants_indices.pl sample.hapcut2.fragment_matrix_file`

The corresponding sample output files are provided in the subdirectory
`sample_output`.

