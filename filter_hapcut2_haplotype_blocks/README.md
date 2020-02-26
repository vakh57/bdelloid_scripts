# Description:
The `filter_hapcut2_haplotype_blocks.pl` script is meant to be used for
filtering of phased haplotype blocks assembled with HapCUT2. The script uses a
list of conflicting variant pairs produced with the
`get_conflicting_variants_indices.pl` script.

The `filter_hapcut2_haplotype_blocks.pl` script reads in two files: 
1. a haplotype file produced with [HapCUT2](https://github.com/vibansal/HapCUT2).
2. a file containing indices of conflicting variant pairs produced with the
`get_conflicting_variants_indices.pl` script.

This list is supposed to be obtained based on the same data as the input
HapCUT2 haplotype file.  i.e. to obtain a list of conflicting variant pairs,
run `get_conflicting_variants_indices.pl` on the same fragment file that was
used as input to HapCUT2 to assemble haplotypes. 
Conflicting variant pairs are defined as pairs of variants for which more than
two haplotypes are found in aligned reads from a single individual. For
details, see README.md provided with the `get_conflicting_variants_indices.pl`
script.

The script filters phased haplotype blocks assembled with HapCUT2 based on: 
* presence of variants involved in conflicting variant pairs. 
* presence of "problematic" variants - these are variants with values of switch
  and/or mismatch quality computed by HapCUT2 below a certain threshold. These
qualities correspond to phred-scaled probabilities of error (see HapCUT2 
[documentation](https://github.com/vibansal/HapCUT2)). The current version of 
the script uses a threshold value of `100`. 
This behaviour can be modified by changing the value assigned to the variable
`$EXP_QUAL_VALUE` inside the script.
Pruned variants are also regarded as problematic. 

The script generates two sets of filtered haplotype blocks:

1. Set 1 includes original phased haplotype blocks with blocks harboring any
phased variant involved in any conflicting pair filtered out. No further
filtering based on switch/mismatch qualities is applied to this set.

2. Set 2 is produced from the Set 1 (that is from the blocks left after removal
of blocks with conflicting variants) based on switch/mismatch qualities. For
each block, the script determines the number of problematic variants (those
with switch/mismatch quality below a certain threshold).  Blocks devoid of
problematic variants are left intact and included in the Set 2. Blocks carrying
more than one problematic variant are discarded. Blocks carrying exactly one
problematic variant are split at the corresponding position. The two chunks
resulting from the split are checked for the number of remaining phased
variants. If a chunk has at least two phased variants, it is retained and
included in the Set 2 as a separate haplotype block.

The threshold value of 100 is used in the current script version. This
threshold can be modified by changing the value assigned to the variable
`$EXP_QUAL_VALUE` inside the script.

# Requirements:

* This script requires Perl 5 (tested with Perl v5.10.1 under Scientific Linux
  release 6.8, Carbon).
* [YAML](https://metacpan.org/pod/YAML) Perl module
* Input haplotype file produced with
  [HapCUT2](https://github.com/vibansal/HapCUT2) in the HapCUT2 format.
* Input file containing indices of conflicting variant pairs produced with the
  `get_conflicting_variants_indices.pl` script from the same fragment file that
was used to assemble haplotypes with HapCUT2. The format of this file is
described in `README.md` provided with the `get_conflicting_variants_indices.pl`
script.

# Usage:

`perl filter_hapcut2_haplotype_blocks.pl input_hapcut2_haplotype_file
input_conflicting_indices_file`

where:
* `input_hapcut2_haplotype_file` is the name of the input haplotype file
  produced with [HapCUT2](https://github.com/vibansal/HapCUT2).
* `input_conflicting_indices_file` is the name of the input file containing
  indices of conflicting variant pairs produced with the
`get_conflicting_variants_indices.pl` script.


# Output files:

The script produces two output files:

* output file with the extension `.conflicts.filtered.txt` contains filtered
  Set 1. This set includes those original intact haplotype blocks that were
free of variants involved in conflicting variant pairs.

* output file with the extension `.conflicts.qual.filtered.txt` contains
  filtered Set 2.  This set includes:
  * original intact haplotype blocks that were free both of variants involved
in conflicts and problematic variants 
  * chunks of original haplotype blocks
resulting from the split of an original block at the position of a problematic
variant

## Format of output files:
Output files are in the same format as an input file with phased haplotype
blocks, that is, in the HapCUT2 haplotype format.  The only difference being
that those haplotype blocks corresponding to chunks resulting from the split of
an original HapCUT2 block, have '0' value in the header field describing the
number of fragments in the block.


# Sample input data:

A sample input haplotype file (`sample.hapcut2.haplotype_file`) and a sample
input file containing indices of conflicting variant pairs
(`sample.conflicting_indices_file`) are included along with the script.

The command to run the script with the provided sample files:
`perl filter_hapcut2_haplotype_blocks.pl sample.hapcut2.haplotype_file
sample.conflicting_indices_file`

The corresponding sample output files are provided in the subdirectory
`sample_output`.

