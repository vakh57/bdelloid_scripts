# Description:

The `find_haplotypic.counterparts.pl` script takes as input reconstructed 
sequences of haplotypes of several individuals from the same locus and 
identifies cases where both haplotypes of the same individual have 
reciprocal closest counterparts in other individuals.

The script is supposed to work with haplotypic sequences reconstructed by 
'inserting' phased SNPs into the corresponding sequence from the reference 
genome. All haplotypes representing the same genomic locus must be  
reconstructed using the same reference sequence and only SNP data (no INDELs)
are assumed to be used for reconstruction. Therefore, the haplotypes obtained 
in this way are effectively 'aligned' and the input sequences are expected 
to be of exactly the same length. 
The script is not supposed to work with haplotypes reconstructed using INDEL
variants. Although, one may try to first align such haplotypes and then use
them as the input data to the script. In any case, gap columns (if present)
will be skipped.

For each of the two haplotypes within each individual present in the input
file, the script computes the nucleotide distance (proportion of nucleotide 
differences) to the other haplotype within the same individual and to each
haplotype in all other individuals. For each haplotype, the script identifies 
the closest haplotypic counterpart in other individuals. 
To test the robustness of this matching, the script compares the number of 
nucleotide differences between the haplotype and its closest (N1) and second 
closest (N2) counterpart. The haplotype is defined as having an unambiguous 
closest counterpart if the difference between N2 and N1 is 3 SNPs or more (N2-N1≥3).
This threshold can be modified by changing the value assigned to the variable
`$MIN_SNP_DIST` inside the script.
Next, for each individual present in the input file, the script checks whether
both its haplotypes (H1 and H2) have unambiguous closest counterparts (H1´ and H2´).
Only such cases are processed further.
Then, the script retains only reciprocal best matches. 
That is, it is required that H1 and H2 are also identified as the closest counterparts 
of the haplotypes H1´ and H2´ respectively.
Such recirpocal groupings of haplotypes are further classified as 'congruent' 
(if reciprocal closest counterparts H1´ and H2´ are found in the same individual)
or 'incongruent' (if reciprocal closest counterparts H1´ and H2´ are found 
in different individuals).


# Requirements:

* This script requires Perl 5 (tested with Perl v5.10.1 under Scientific Linux
  release 6.8, Carbon).
* [File::Basename](https://perldoc.perl.org/File/Basename.html) Perl module.
* Input file containing data on haplotypes from several individuals in the FASTA format. 
  Each individual is expected to be represented by two haplotypes. The two haplotypes
  of the same individual have to be designated as '.hap1' and '.hap2'
  (e.g. Ind1.hap1 and Ind1.hap2). 
  All haplotypes contained in the input FASTA file are expected to be obtained by 
  'inserting' phased SNPs in the same reference sequence. That is, all haplotypes
  have to be in the same coordinate system ('aligned'). Alternatively, haplotypes
  reconstructed using both SNPs and INDELs have to be aligned before running the 
  script. Gap columns will be skipped, if present.

# Usage:

`perl find_haplotypic.counterparts.pl input_haplotypes.fasta`

where:
* `input_haplotypes.fasta` is the name of the input FASTA file 
 containing reconstructed sequences of haplotypes of several 
 individuals from the same genomic locus.

# Output files:

The script produces three output files (the second and the third output 
files contain subsets of the data contained in the first output file):

* output file with the extension `.recipr.counterparts.txt` contains data
  on reciprocal closest counterparts identified for the haplotypes present
  in the input FASTA file. Only those cases where reciprocal closest 
  counterparts have been identified for both haplotypes of an individual are 
  considered. Both congruent and incongruent haplotype groupings are included.

* output file with the extension `.recipr.counterparts.diff_ind.txt` contains
  only incongruent haplotype groupings. That is, only such cases where both 
  haplotypes of an individual have reciprocal closest counterparts and these
  counterparts are found in different individuals are listed.

* output file with the extension `.recipr.counterparts.same_ind.txt` contains
  only congruent haplotype groupings. That is, only such cases where both 
  haplotypes of an individual have reciprocal closest counterparts and these
  counterparts are found in the same individual are listed.

If an output file is empty, this means that groupings of the corresponding
type have not been found.


## Format of output files:

Output files are tab-delimited.
Each line of the output file contains information on reciprocal haplotypic
counterparts for the two haplotypes of the same individual.

Field 1 is the name of the input FASTA file.
Field 2 is the ID of individual.
Field 3 is the pattern of groupings found for the two haplotypes of this
individual.
Field 4 is the ID of the first haplotype of this individual (H1).
Field 5 is the ID of the reciprocal haplotypic counterpart (H1´) identified for H1.
Field 6 is the nucleotide distance (proportion of nucleotide differences) 
between H1 and H1´.
Field 7 is the absolute nucleotide distance (absolute number of nucleotide 
differences) between H1 and H1´.
Field 8 is the ID of the second haplotype of this individual (H2).
Field 9 is the ID of the reciprocal haplotypic counterpart (H2´) identified for H2.
Field 10 is the nucleotide distance (proportion of nucleotide differences) 
between H2 and H2´.
Field 11 is the absolute nucleotide distance (absolute number of nucleotide 
differences) between H2 and H2´.
Field 12 is the total number of analyzed sites (free of gaps and ambiguous symbols).
Field 13 is the total number of sites in the input 'alignment'. 
Field 14 is the duplicate of the field 2 (ID of individual).
Field 15 is the nucleotide distance (proportion of nucleotide differences)
between the two haplotypes of this individual (H1 and H2). 
This is intraindividual distance.
Field 16 is the absolute nucleotide distance (absolute number of nucleotide
differences) between H1 and H2.

# Sample input data:

A sample input file (`sample.haps.fasta`) 
containing sequences of haplotypes representing the same locus 
in 8 individuals is included along with the script.

The command to run the script with the sample input file:
`perl find_haplotypic.counterparts.pl sample.haps.fasta`

The corresponding sample output files are provided in the subdirectory
`sample_output`.

