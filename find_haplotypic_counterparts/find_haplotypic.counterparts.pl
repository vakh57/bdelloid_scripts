#!/usr/bin/env perl

# This file is part of the software that was used to process sequencing data 
# for A. vaga individuals L1-L11 analyzed in the manuscript 
# https://www.biorxiv.org/content/10.1101/489393v1 
#
# Copyright (c) 2018-2020, by Olga Vakhrusheva <vakh57@gmail.com>
#
# The script is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# The script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with the software; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.

use strict;
use List::Util qw(min max sum);
use Cwd;
use File::Basename;

my ($fasta_file) = @ARGV;

unless (defined($fasta_file))
{
    die "\nInput file is not specified. Stopped.";
}

print "\n";
print "#######################################################################################################\n";
print "Run find_haplotypic_counterparts.pl\n";
print "#######################################################################################################\n";

my $cur = cwd();
# print "$cur\n";

my $file_id = basename($fasta_file);
# print "file_id: $file_id\n";
# print "\n";

my $file_path = dirname($fasta_file);
# print "\n";
# print "file_path: $file_path\n";

(my $file_base = $file_id) =~ s/\.fasta$//;

# print "file_base: $file_base\n";

print "\n";
print "This script takes as input reconstructed sequences of haplotypes of several individuals from the same\nlocus and identifies cases where both haplotypes of the same individual have reciprocal closest counterparts.\n";
print "\n";
print "For each of the two haplotypes within each individual present in the input file, the script computes\nthe nucleotide distance (proportion of nucleotide differences) to the other haplotype within the same\nindividual and to each haplotype in all other individuals.\nFor each haplotype, the script identifies the closest haplotypic counterpart in other individuals.\n";
print "To test the robustness of this matching, the script compares the number of nucleotide differences\nbetween the haplotype and its closest (N1) and second closest (N2) counterpart.\n";
print "The haplotype is defined as having an unambiguous closest counterpart if the difference between\nN2 and N1 is 3 SNPs or more (N2-N1≥3).\n";
print "This threshold can be modified by changing the value assigned to the variable \$MIN_SNP_DIST\ninside the script.\n";
print "Next, for each individual present in the input file, the script checks whether both its haplotypes\n(H1 and H2) have unambiguous closest counterparts (H1´ and H2´).\n";
print "Only such cases are processed further.\n";
print "Then, the script retains only reciprocal best matches.\n";
print "That is, it is required that H1 and H2 are also identified as the closest counterparts of the\nhaplotypes H1´ and H2´ respectively.\n";
print "Such recirpocal groupings of haplotypes are further classified as 'congruent' (if reciprocal\nclosest counterparts H1´ and H2´ are found in the same individual) or 'incongruent'\n(if reciprocal closest counterparts H1´ and H2´ are found in different individuals).\n";
print "\n";
print "For further information, see README.md.\n";


#######################################################################################################
### Specify output file names

### File to report data on reciprocal haplotypic counterparts 
### Only those cases when both haplotypes of an individual have unambigous closest counterparts are reported
my $CLNEIGHBOR_FILE = "$file_base.recipr.counterparts.txt";
open (CLFILE,">$CLNEIGHBOR_FILE") or die $!;

### File to report data on reciprocal haplotypic counterparts found in different individuals
### Only those cases when both haplotypes of an individual have unambigous closest counterparts are reported
### Among those, the cases where the counterparts of the first and the second haplotype are found in different individuals are selected
my $CLNEIGHBOR_FILE_DIND = "$file_base.recipr.counterparts.diff_ind.txt";
open (CLFILE_DIFF_IND,">$CLNEIGHBOR_FILE_DIND") or die $!;

### File to report data on reciprocal haplotypic counterparts found in the same individual
### Only those cases when both haplotypes of an individual have unambigous closest counterparts are reported
### Among those, the cases where the counterparts of the first and the second haplotype are found in the same individual are selected
my $CLNEIGHBOR_FILE_SIND = "$file_base.recipr.counterparts.same_ind.txt";
open (CLFILE_SAME_IND,">$CLNEIGHBOR_FILE_SIND") or die $!;

my $sepstring = '#######################################################################################################';

print "\n";
print "#######################################################################################################\n";
print "\n";
print "Started\n";
system("date");
print "\n";
print "#######################################################################################################\n";

print "\n";
print "#######################################################################################################\n";
print "#######################################################################################################\n";
print "Read the input FASTA file:\n";
print "\n";
print "Input FASTA file:\n$fasta_file\n";
print "\n";

print "Number of entries in the input FASTA file (presumed haplotypes):\n";
system("grep -c '^>' $fasta_file");


#######################################################################################################
### Proceed to reading the input FASTA file presumed to contain sequences of haplotypes from several individuals
open (MF,"$fasta_file") or die $!;

my $seq_id;

### Reference to a hash used to store sequences of haplotypes
my $seq_hash;

while (<MF>)
{
    chomp();
### Process FASTA description lines
    if (/^>/)
    {
	$seq_id = $_;
	$seq_id =~ s/^>//;
    }
### Process FASTA lines containing sequence data
    else
    {

### Check that the sequence ID is set
	unless ($seq_id)
	{
	    die "\nInput file doesn't conform to the FASTA format. Stopped.";
	}

	$seq_hash->{$seq_id}.=$_;
    }

}

close(MF);

print "#######################################################################################################\n";

### Retrieve indices of sites free of gaps and ambiguous symbols (other than A, T, G or C)
my ($locus_unambiguous_sites_ref, $tot_unambiguous_sites, $tot_sites) = find_ungapped_sites($seq_hash);

### Warn if the number of ungapped columns < $MIN_SITES
my $MIN_SITES = 100;

my $total_unambiguous_sites_ref = scalar(@{$locus_unambiguous_sites_ref});

if ($total_unambiguous_sites_ref < $MIN_SITES)
{
    warn "\nWARNING: Too few sites without gaps are available for the analysis in $file_id: $total_unambiguous_sites_ref\n";
}

print "\nTotal number of sites available for the analysis: $total_unambiguous_sites_ref\n";
print "#######################################################################################################\n";

### Retrieve IDs of individuals from IDs of haplotypes present in the input FASTA file
### Two haplotypes of the same individual (IND) must be designated as IND.hap1 and IND.hap2
### Store IDs of individuals in the hash

my $sample_ids_hash;

foreach my $seqkey (sort {$a cmp $b} keys(%{$seq_hash}))
{

    if ($seqkey =~ /^(\S+)\.hap[12]$/)
    {
	# print "Haplotype ID $seqkey -> Individual ID $1\n"; 
	$sample_ids_hash->{$1} = 1;
    }
    else
    {
	die "\nHaplotype ID does not conform to the expected format: $seqkey\nTwo haplotypes of the same individual (IND) must be designated as IND.hap1 and IND.hap2.\nStopped.\n";
    }

}

my @predefined_sample_ids = (sort {$a cmp $b} keys(%{$sample_ids_hash}));
my $npredefined_samples = scalar(@predefined_sample_ids);

my $nexp_samples_to_compare = $npredefined_samples-1;

print "\n";
print "IDs of individuals whose haplotypes are present in the input FASTA file:\n";
print "@predefined_sample_ids\n";
print "$sepstring\n";


#######################################################################################################
#######################################################################################################
### Iterate over individuals and identify closest counterparts for both haplotypes of each individual
### A haplotype is then defined as having an unambiguous closest counterpart, if the difference in distances from this haplotype to its closest and second closest counterparts >= $MIN_SNP_DIST

#######################################################################################################
### Compute intraindividual haplotypic distances (distances between the two haplotypes of the same individual)

my $WITHIN_SAMPLE_HAPDISTANCES;
my $WITHIN_SAMPLE_HAP_ABSDISTANCES;

### Iterate over individuals
foreach my $present_sampleid (@predefined_sample_ids)
{

    my $sample_al1 = "$present_sampleid.hap1";
    my $sample_al2 = "$present_sampleid.hap2";

    my ($diverged_fraction_new,$ndiverged_sites_new,$nlookedup_sites_new) = get_pairwise_distance($seq_hash,$locus_unambiguous_sites_ref,$sample_al1,$sample_al2);

### Store the distance expressed as a proportion of nucleotide differences
    $WITHIN_SAMPLE_HAPDISTANCES->{$present_sampleid} = $diverged_fraction_new;

### Store the absolute distance expressed as an absolute number of nucleotide differences
    $WITHIN_SAMPLE_HAP_ABSDISTANCES->{$present_sampleid} = $ndiverged_sites_new;

    # print "$present_sampleid: $diverged_fraction_new $ndiverged_sites_new\n";
}


#######################################################################################################
### Compute distances between the two haplotypes of each individual and haplotypes in the rest of individuals

# print "\n";
# print "Compute distances between the two haplotype of each individual and haplotypes in the rest of individuals\n";

my $MIN_SNP_DIST = 3;

my $ALL_SAMPLE_PAIRS_IDS;
my $ALL_SAMPLE_PAIRS_ABS_IDS;

### Reference to a hash used to index cases when at least one haplotype has an unambiguous closest counterpart
my $CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX; 

for (my $ind0=0;$ind0<@predefined_sample_ids;$ind0++)
{
    my $given_sample1 = $predefined_sample_ids[$ind0];

### IDs of the two haplotypes of the first individual in the pair
    my @sample1_alleles = map { my $alid="$given_sample1.$_"; $alid; } qw (hap1 hap2); 

    my $samp1_al_first = "$given_sample1.hap1";
    my $samp1_al_second = "$given_sample1.hap2";

### References to hashes used to store distances between the two haplotypes of the current individual and the rest of haplotypes from other individuals
    my $given_sample_ids_hash;
    my $given_sample_abs_ids_hash;

    for (my $ind1=0;$ind1<@predefined_sample_ids;$ind1++)
    {

	my $given_sample2 = $predefined_sample_ids[$ind1];

### IDs of the two haplotypes of the second individual in the pair
	my @sample2_alleles = map { my $alid="$given_sample2.$_"; $alid; } qw (hap1 hap2); 

	my $samp2_al_first = "$given_sample2.hap1";
	my $samp2_al_second = "$given_sample2.hap2";

	foreach my $samp1_al (@sample1_alleles)
	{
	    my @samp1_al_abs_dists;

	    foreach my $samp2_al (@sample2_alleles)
	    {

		my ($diverged_fraction_pair,$ndiverged_sites_pair,$nlookedup_sites_pair) = get_pairwise_distance($seq_hash,$locus_unambiguous_sites_ref,$samp1_al,$samp2_al);

		push @samp1_al_abs_dists,$ndiverged_sites_pair;

		$given_sample_ids_hash->{$samp1_al}->{$samp2_al} = $diverged_fraction_pair;

		$given_sample_abs_ids_hash->{$samp1_al}->{$samp2_al} = $ndiverged_sites_pair;

		$ALL_SAMPLE_PAIRS_IDS->{$given_sample1}->{$given_sample2}->{$samp1_al}->{$samp2_al} = $diverged_fraction_pair;
		$ALL_SAMPLE_PAIRS_ABS_IDS->{$given_sample1}->{$given_sample2}->{$samp1_al}->{$samp2_al} = $ndiverged_sites_pair;

	    }

	}

    }



#######################################################################################################
### Look at the distances from the two haplotypes of the current individual ($given_sample1) to haplotypes from the rest of individuals

### Array containing IDs of the individuals other than the current one ($given_sample1)
    my @samples_to_compare = grep { $_ ne $given_sample1; } @predefined_sample_ids;
    my $dim_samples_to_compare = scalar(@samples_to_compare);

    unless ($dim_samples_to_compare == $nexp_samples_to_compare)
    {
	die "\nUnexpected number of individuals to look at $given_sample1: $dim_samples_to_compare vs $nexp_samples_to_compare. Stopped.\n";
    }

### Array containing IDs of first haplotypes (.hap1) carried by the rest of individuals
    my @alids1_to_compare = map { my $alid1="$_.hap1"; $alid1; } @samples_to_compare;

### Array containing IDs of second haplotypes (.hap2) carried by the rest of individuals
    my @alids2_to_compare = map { my $alid2="$_.hap2"; $alid2; } @samples_to_compare;

    my @alidsboth_to_compare = (@alids1_to_compare,@alids2_to_compare);

### Check whether all distances that were supposed to be computed were computed and stored in the hash
### These are distances from both haplotypes of this individual ($given_sample1) to the rest of haplotypes in other individuals
    my @alids_check_array = map {
	if ((exists ($given_sample_ids_hash->{$samp1_al_first}->{$_})) && (exists ($given_sample_ids_hash->{$samp1_al_second}->{$_}) ))
	{ 1; }
	else 
	{ 0; }

    } @alidsboth_to_compare;

    my $sum_alids_check_array = sum(@alids_check_array);
    
    # print "sum_alids_check_array $sum_alids_check_array\n";
    # print "\n";

    unless ($sum_alids_check_array == ($nexp_samples_to_compare*2))
    {
	die "\nSome of the data on interhaplotype distances for $given_sample1 are missing: only $sum_alids_check_array distances are set. Stopped.\n";

    }


### Sort haplotypes from other individuals based on their distances to the first haplotype of the current individual ($given_sample1)
    my @samp1_al_first_sorted_haplotypes = sort { ($given_sample_ids_hash->{$samp1_al_first}->{$a} <=> $given_sample_ids_hash->{$samp1_al_first}->{$b}) || ($a cmp $b) } @alidsboth_to_compare;
    my @samp1_al_first_sorted_haplotypes_dist = map { $given_sample_ids_hash->{$samp1_al_first}->{$_}; } @samp1_al_first_sorted_haplotypes;

    # print "\n";
    # print "First allele:\n";
    # print "$samp1_al_first -> samp1_al_first_sorted_haplotypes:\n";
    # print "@samp1_al_first_sorted_haplotypes\n";
    # print "samp1_al_first_sorted_haplotypes_dist:\n";  
    # print "@samp1_al_first_sorted_haplotypes_dist\n";

### Closest counterpart of the first haplotype and the corresponding distance 
    my $samp1_al_first_neighbor1 = $samp1_al_first_sorted_haplotypes[0];
    my $samp1_al_first_neighbor1_dist = $given_sample_ids_hash->{$samp1_al_first}->{$samp1_al_first_neighbor1};
    my $samp1_al_first_neighbor1_abs_dist = $given_sample_abs_ids_hash->{$samp1_al_first}->{$samp1_al_first_neighbor1};

### Determine which individual carries the closest counterpart of the first haplotype ($samp1_al_first_neighbor1)
    my $samp1_al_first_neighbor1_sample;

    if ($samp1_al_first_neighbor1 =~ /^(\S+)\.hap[12]$/)
    {
	$samp1_al_first_neighbor1_sample = $1;
    }
    else
    {
	die "\nID of the putative haplotypic counterpart for $samp1_al_first does not conform to the expected format: $samp1_al_first_neighbor1. Stopped.\n";
    }

### Second closest counterpart of the first haplotype and the corresponding distance
    my $samp1_al_first_neighbor2 = $samp1_al_first_sorted_haplotypes[1];
    my $samp1_al_first_neighbor2_abs_dist = $given_sample_abs_ids_hash->{$samp1_al_first}->{$samp1_al_first_neighbor2};

    # print "samp1_al_first_neighbor1: $samp1_al_first_neighbor1 -> $samp1_al_first_neighbor1_dist ($samp1_al_first_neighbor1_abs_dist)\n";
    # print "samp1_al_first_neighbor1_sample: $samp1_al_first_neighbor1_sample\n";
    # print "\n";
    # print "samp1_al_first_neighbor1: $samp1_al_first_neighbor1 -> $samp1_al_first_neighbor1_abs_dist\n";
    # print "samp1_al_first_neighbor2: $samp1_al_first_neighbor2 -> $samp1_al_first_neighbor2_abs_dist\n";

### Sort haplotypes from other individuals based on their distances to the second haplotype of the current individual ($given_sample1)
    my @samp1_al_second_sorted_haplotypes = sort { ($given_sample_ids_hash->{$samp1_al_second}->{$a} <=> $given_sample_ids_hash->{$samp1_al_second}->{$b}) || ($a cmp $b) } @alidsboth_to_compare;
    my @samp1_al_second_sorted_haplotypes_dist = map { $given_sample_ids_hash->{$samp1_al_second}->{$_}; } @samp1_al_second_sorted_haplotypes;

    # print "\n";
    # print "Second allele:\n";
    # print "$samp1_al_second -> samp1_al_second_sorted_haplotypes:\n";
    # print "@samp1_al_second_sorted_haplotypes\n";
    # print "samp1_al_second_sorted_haplotypes_dist:\n";  
    # print "@samp1_al_second_sorted_haplotypes_dist\n";

### Closest counterpart of the second haplotype and the corresponding distance
    my $samp1_al_second_neighbor1 = $samp1_al_second_sorted_haplotypes[0];
    my $samp1_al_second_neighbor1_dist = $given_sample_ids_hash->{$samp1_al_second}->{$samp1_al_second_neighbor1};
    my $samp1_al_second_neighbor1_abs_dist = $given_sample_abs_ids_hash->{$samp1_al_second}->{$samp1_al_second_neighbor1};

### Determine which individual carries the closest counterpart of the second haplotype ($samp1_al_second_neighbor1)
    my $samp1_al_second_neighbor1_sample;

    if ($samp1_al_second_neighbor1 =~ /^(\S+)\.hap[12]$/)
    {
	$samp1_al_second_neighbor1_sample = $1;
    }
    else
    {
	die "\nID of the putative haplotypic counterpart for $samp1_al_second does not conform to the expected format: $samp1_al_second_neighbor1. Stopped.\n";

    }

### Second closest counterpart of the second haplotype and the corresponding distance
    my $samp1_al_second_neighbor2 = $samp1_al_second_sorted_haplotypes[1];
    my $samp1_al_second_neighbor2_abs_dist = $given_sample_abs_ids_hash->{$samp1_al_second}->{$samp1_al_second_neighbor2};

    # print "samp1_al_second_neighbor1: $samp1_al_second_neighbor1 -> $samp1_al_second_neighbor1_dist ($samp1_al_second_neighbor1_abs_dist)\n";
    # print "samp1_al_second_neighbor1_sample: $samp1_al_second_neighbor1_sample\n";
    # print "\n";
    # print "samp1_al_second_neighbor1: $samp1_al_second_neighbor1 -> $samp1_al_second_neighbor1_abs_dist\n";
    # print "samp1_al_second_neighbor2: $samp1_al_second_neighbor2 -> $samp1_al_second_neighbor2_abs_dist\n";

### Check whether the difference in distances from the first haplotype to its second and first closest counterparts >= $MIN_SNP_DIST
### If this condition is met, remember that the first haplotype has an unambiguous closest counterpart
    if (($samp1_al_first_neighbor2_abs_dist-$samp1_al_first_neighbor1_abs_dist) >= $MIN_SNP_DIST) 
    {
	$CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$given_sample1}->{$samp1_al_first} = $samp1_al_first_neighbor1;
    }

### Check whether the difference in distances from the second haplotype to its second and first closest counterparts >= $MIN_SNP_DIST
### If this condition is met, remember that the second haplotype has an unambiguous closest counterpart
    if (($samp1_al_second_neighbor2_abs_dist-$samp1_al_second_neighbor1_abs_dist) >= $MIN_SNP_DIST)
    {
	$CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$given_sample1}->{$samp1_al_second} = $samp1_al_second_neighbor1;
    }

}

#######################################################################################################
#######################################################################################################
### Check whether the identified haplotypic counterparts are reciprocal best matches\n
print "\n";
print "\n";
print "$sepstring\n";
print "$sepstring\n";
print "Look for cases when unambiguous closest counterparts of the two haplotypes correspond to reciprocal best matches\n";


### Iterate over all individuals for which at least one haplotype has an unambiguous closest counterpart
CHECKSAMPLES: foreach my $NSAMPLE (sort { $a cmp $b } keys(%{$CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX}))
{
    my $NSAMPLE_AL1 = "$NSAMPLE.hap1";
    my $NSAMPLE_AL2 = "$NSAMPLE.hap2";

### Retrieve the distance between the two haplotypes of the current individual ($NSAMPLE)
### This is an intraindividual distance (between the two haplotypes of the same individual)
    my $NSAMPLE_INTRA_HAP_DIST = $WITHIN_SAMPLE_HAPDISTANCES->{$NSAMPLE};
    my $NSAMPLE_INTRA_HAP_DIST_ABS = $WITHIN_SAMPLE_HAP_ABSDISTANCES->{$NSAMPLE};

### Consider further only those cases when unambiguous closest counterparts in other individuals were identified for both haplotypes of the current individual 

    if ((exists($CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$NSAMPLE}->{$NSAMPLE_AL1})) && (exists($CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$NSAMPLE}->{$NSAMPLE_AL2})))
    {

	print "$sepstring\n";
	print "\nBoth haplotypes of the individual $NSAMPLE have unambiguous closest counterparts.\n";

### Closest counterpart of the first haplotype of the individual $NSAMPLE
	my $NSAMPLE_AL1_CLN = $CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$NSAMPLE}->{$NSAMPLE_AL1};

### Determine which individual carries the closest counterpart ($NSAMPLE_AL1_CLN) of the first haplotype 
	my $NSAMPLE_AL1_CLN_SAMPLE;

	if ($NSAMPLE_AL1_CLN =~ /^(\S+)\.hap[12]$/)
	{
	    $NSAMPLE_AL1_CLN_SAMPLE = $1;
	}
	else
	{
	    die "\nID of the putative haplotypic counterpart for $NSAMPLE_AL1 does not conform to the expected format: $NSAMPLE_AL1_CLN. Stopped.\n";
	}

### Retrieve the distance from the first haplotype of the considered individual to its closest counterpart in another individual
	my $NSAMPLE_AL1_CLN_DIST = $ALL_SAMPLE_PAIRS_IDS->{$NSAMPLE}->{$NSAMPLE_AL1_CLN_SAMPLE}->{$NSAMPLE_AL1}->{$NSAMPLE_AL1_CLN};
	my $NSAMPLE_AL1_CLN_ABS_DIST = $ALL_SAMPLE_PAIRS_ABS_IDS->{$NSAMPLE}->{$NSAMPLE_AL1_CLN_SAMPLE}->{$NSAMPLE_AL1}->{$NSAMPLE_AL1_CLN};

### Retrieve the distance between the two haplotypes within the individal ($NSAMPLE_AL1_CLN_SAMPLE) harbouring the closest counterpart ($NSAMPLE_AL1_CLN) of the first haplotype ($NSAMPLE_AL1)
### This is an intraindividual distance (between the two haplotypes of the same individual)
	my $NSAMPLE_AL1_CLN_SAMPLE_INTRA_HAP_DIST = $WITHIN_SAMPLE_HAPDISTANCES->{$NSAMPLE_AL1_CLN_SAMPLE};

### Closest counterpart of the second haplotype of the individual $NSAMPLE
	my $NSAMPLE_AL2_CLN = $CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$NSAMPLE}->{$NSAMPLE_AL2};

### Determine which individual carries the closest counterpart ($NSAMPLE_AL2_CLN) of the second haplotype
	my $NSAMPLE_AL2_CLN_SAMPLE;

	if ($NSAMPLE_AL2_CLN =~ /^(\S+)\.hap[12]$/)
	{
	    $NSAMPLE_AL2_CLN_SAMPLE = $1;
	}
	else
	{
	    die "\nID of the putative haplotypic counterpart for $NSAMPLE_AL2 does not conform to the expected format: $NSAMPLE_AL2_CLN. Stopped.\n";
	}

### Retrieve the distance from the second haplotype of the considered individual to its closest counterpart in another individual
	my $NSAMPLE_AL2_CLN_DIST = $ALL_SAMPLE_PAIRS_IDS->{$NSAMPLE}->{$NSAMPLE_AL2_CLN_SAMPLE}->{$NSAMPLE_AL2}->{$NSAMPLE_AL2_CLN};
	my $NSAMPLE_AL2_CLN_ABS_DIST = $ALL_SAMPLE_PAIRS_ABS_IDS->{$NSAMPLE}->{$NSAMPLE_AL2_CLN_SAMPLE}->{$NSAMPLE_AL2}->{$NSAMPLE_AL2_CLN};

### Retrieve the distance between the two haplotypes within the individal ($NSAMPLE_AL2_CLN_SAMPLE) harbouring the closest counterpart ($NSAMPLE_AL2_CLN) of the second haplotype ($NSAMPLE_AL2)
### This is an intraindividual distance (between the two haplotypes of the same individual)
	my $NSAMPLE_AL2_CLN_SAMPLE_INTRA_HAP_DIST = $WITHIN_SAMPLE_HAPDISTANCES->{$NSAMPLE_AL2_CLN_SAMPLE};

	# print "Haplotype 1: $NSAMPLE_AL1 -> $NSAMPLE_AL1_CLN (distance $NSAMPLE_AL1_CLN_DIST) ($NSAMPLE_AL1_CLN_SAMPLE $NSAMPLE_AL1_CLN_SAMPLE_INTRA_HAP_DIST)\n";
	# print "Haplotype 2: $NSAMPLE_AL2 -> $NSAMPLE_AL2_CLN (distance $NSAMPLE_AL2_CLN_DIST) ($NSAMPLE_AL2_CLN_SAMPLE $NSAMPLE_AL2_CLN_SAMPLE_INTRA_HAP_DIST)\n";

	print "Haplotype 1: $NSAMPLE_AL1 -> $NSAMPLE_AL1_CLN\n";
	print "Haplotype 2: $NSAMPLE_AL2 -> $NSAMPLE_AL2_CLN\n";

	my @SORTED_UNAMBIG_NEIGHBORS = sort { $a cmp $b } ($NSAMPLE_AL1_CLN_SAMPLE, $NSAMPLE_AL2_CLN_SAMPLE);
	my $SORTED_UNAMBIG_NEIGHBORS_STRING = join('.vs.', @SORTED_UNAMBIG_NEIGHBORS);

	my $LINE_CNEIGHBORS = "$file_id\t$NSAMPLE\t$SORTED_UNAMBIG_NEIGHBORS_STRING\t$NSAMPLE_AL1\t$NSAMPLE_AL1_CLN\t$NSAMPLE_AL1_CLN_DIST\t$NSAMPLE_AL1_CLN_ABS_DIST\t$NSAMPLE_AL2\t$NSAMPLE_AL2_CLN\t$NSAMPLE_AL2_CLN_DIST\t$NSAMPLE_AL2_CLN_ABS_DIST\t$tot_unambiguous_sites\t$tot_sites\t$NSAMPLE\t$NSAMPLE_INTRA_HAP_DIST\t$NSAMPLE_INTRA_HAP_DIST_ABS";

	if (($NSAMPLE_AL1_CLN_SAMPLE ne $NSAMPLE) && ($NSAMPLE_AL2_CLN_SAMPLE ne $NSAMPLE))
	{
	    my $DIST_WEIRD = 0;

### Check that the distance between the two haplotypes within the considered individals is larger than distances from both haplotypes to the closest neighbors in other individuals 
	    unless (($NSAMPLE_AL1_CLN_DIST < $NSAMPLE_INTRA_HAP_DIST) && ($NSAMPLE_AL2_CLN_DIST < $NSAMPLE_INTRA_HAP_DIST))
	    {
		$DIST_WEIRD++;
	    }
	    unless (($NSAMPLE_AL1_CLN_DIST < $NSAMPLE_AL1_CLN_SAMPLE_INTRA_HAP_DIST) && ($NSAMPLE_AL2_CLN_DIST < $NSAMPLE_AL2_CLN_SAMPLE_INTRA_HAP_DIST))
	    {
		$DIST_WEIRD++;
	    }

	    if ($DIST_WEIRD > 0)
	    {
		# print "\n";
		# print "$NSAMPLE did not pass check on distances between haplotypes, proceed to next the sample.\n";
		next CHECKSAMPLES;
	    }
	    else
	    {
		# print "\n";
                # print "$NSAMPLE passed check on distances between haplotypes.\n";
		;
### Check whether the counterparts for the two haplotypes are found in different individuals
		# if ($NSAMPLE_AL1_CLN_SAMPLE ne $NSAMPLE_AL2_CLN_SAMPLE)
		# {
		#     ;
		# }

	    }

### Check whether the identified counterparts are reciprocal best matches
	    if ((exists($CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$NSAMPLE_AL1_CLN_SAMPLE}->{$NSAMPLE_AL1_CLN})) && (exists($CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$NSAMPLE_AL2_CLN_SAMPLE}->{$NSAMPLE_AL2_CLN})))
	    {

### Closest counterpart for the closest counterpart of the first haplotype
		my $NSAMPLE_AL1_CLN_CLN = $CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$NSAMPLE_AL1_CLN_SAMPLE}->{$NSAMPLE_AL1_CLN};
		# print "\n";
		# print "Check 1: NSAMPLE_AL1_CLN $NSAMPLE_AL1_CLN -> $NSAMPLE_AL1_CLN_CLN\n";

### Closest counterpart for the closest counterpart of the second haplotype
		my $NSAMPLE_AL2_CLN_CLN = $CLOSEST_SEP_UNAMBIG_NEIGHBORS_INDEX->{$NSAMPLE_AL2_CLN_SAMPLE}->{$NSAMPLE_AL2_CLN};
		# print "Check 2: NSAMPLE_AL2_CLN $NSAMPLE_AL2_CLN -> $NSAMPLE_AL2_CLN_CLN\n";

### Proceed further, only if matches are reciprocal
		if (($NSAMPLE_AL1_CLN_CLN eq $NSAMPLE_AL1) && ($NSAMPLE_AL2_CLN_CLN eq $NSAMPLE_AL2))
		{
		    print "\nThe closest counterparts for haplotypes of the individual $NSAMPLE correspond to reciprocal best matches.\n";

		    print CLFILE "$LINE_CNEIGHBORS\n";

### If counterparts for the two haplotypes are found in different individuals
		    if ($NSAMPLE_AL1_CLN_SAMPLE ne $NSAMPLE_AL2_CLN_SAMPLE)
		    {

			print "\nReciprocal closest counterparts for the two haplotypes of the individual $NSAMPLE are found in different individuals:\n";
			print "$NSAMPLE_AL1_CLN_SAMPLE and $NSAMPLE_AL2_CLN_SAMPLE\n";

			print CLFILE_DIFF_IND "$LINE_CNEIGHBORS\n";

		    }

### If counterparts for the two haplotypes are found in the same individual
		    elsif ($NSAMPLE_AL1_CLN_SAMPLE eq $NSAMPLE_AL2_CLN_SAMPLE)
		    {
			print CLFILE_SAME_IND "$LINE_CNEIGHBORS\n";
		    }

		}

	    } 
	    else
	    {
		print "\nThere is no unambiguous closest counterpart for $NSAMPLE_AL1_CLN or $NSAMPLE_AL2_CLN, or both.\n";
	    }

	}
	else 
	{
	    die "\nCounterparts for haplotypes of $NSAMPLE are found in $NSAMPLE. This should not happen. Stopped.\n";
	}

    }

}


close(CLFILE);
close(CLFILE_DIFF_IND);
close(CLFILE_SAME_IND);


print "#######################################################################################################\n";
print "#######################################################################################################\n";
print "\nFinished\n";
system("date");


### Subroutines

sub find_ungapped_sites
{
    my ($seq_hash_find) = @_;

    my @seqnames = (sort {$a cmp $b} keys(%{$seq_hash_find}));

    print "\nSequence IDs present in the input FASTA file:\n@seqnames\n\n";

    my $first_lineage = $seqnames[0];

### All sequences are assumed to be of the same length
    my $alignment_length = length($seq_hash_find->{$first_lineage});

    my $total_gapped_columns = 0;
    my $total_ungapped_columns = 0;
    my $total_unambiguous_columns = 0;

### Array to store indices of alignment columns (sites) free of gaps
    my @ungapped_columns;

### Array to store indices of alignment columns (sites) free both of gaps and ambiguous characters
    my @unambiguous_columns;

### Iterate over alignment columns (sites)
    for (my $ali=0;$ali<$alignment_length;$ali++)
    {
	# my $alireal = $ali+1;

### Variable used to store the sequence of the current alignment column
	my $aligned_column;

### Iterate over all haplotypes present in the alignment and extract the symbol carried by this haplotype at the current column 
	foreach my $strainal (@seqnames)
	{
	    my $strainal_char = substr($seq_hash_find->{$strainal},$ali,1);

	    $aligned_column.=$strainal_char;
	}

	if ($aligned_column =~ /-/)
	{
	    $total_gapped_columns++;
	}
	else 
	{
	    $total_ungapped_columns++;
	
	    push @ungapped_columns,$ali;
	}

	unless ($aligned_column =~ /[^atgc]/i)
	{
	    $total_unambiguous_columns++;

	    push @unambiguous_columns,$ali;
	}

    }

    unless ($total_unambiguous_columns == $total_ungapped_columns)
    {
    	warn("\nWARNING: Total number of columns free of gaps and ambiguous symbols is not equal to the total number of columns without gaps ($total_unambiguous_columns vs $total_ungapped_columns). Some columns may contain ambiguity characters. This shouldn't happen.\n");	
    }

    print "Total length of the alignment: $alignment_length\n";
    print "\nTotal number of columns with gaps: $total_gapped_columns\n";
    print "Total number of columns without gaps: $total_ungapped_columns\n";
    print "Total number of columns without gaps and ambiguous symbols: $total_unambiguous_columns\n";

    return(\@unambiguous_columns,$total_unambiguous_columns,$alignment_length);

}



sub get_pairwise_distance
{

    my ($seq_hash_local,$fine_sites_hash,$sample1,$sample2) = @_;

### Retrieve the indices of the alignment columns to be used to compute pairwise distance
### These indices are supposed to correspond to columns free of gaps and ambiguity characters
    my @fine_sites_array = @{$fine_sites_hash};
    # print "fine_sites_array: @fine_sites_array\n";

    unless (exists($seq_hash_local->{$sample1}) && (exists($seq_hash_local->{$sample2})))
    {
	die "\nError in the subroutine get_pairwise_distance: one or both sequences with the IDs provided to the subroutine ($sample1 $sample2) are missing from the hash containing the sequence data.\nStopped.\n";
    }

    my $lsample1 = length($seq_hash_local->{$sample1});
    my $lsample2 = length($seq_hash_local->{$sample2});

    unless ($lsample1 == $lsample2)
    {
	die "\nError in the subroutine get_pairwise_distance: sequence lengths (IDs $sample1 and $sample2) do not match: $lsample1 vs $lsample2.\nSequences must be of equal length. Stopped.\n";
    }

    my $nlookedup_sites = 0;
    my $ndiverged_sites = 0;

    foreach my $looksite (sort {$a <=> $b} @fine_sites_array)
    {
	if (($looksite >= $lsample1) || ($looksite < 0))
	{
	    die "\nError in the subroutine get_pairwise_distance: site with the index $looksite is out of the allowed range: sequence length $lsample1.\nStopped.\n";
	}

	$nlookedup_sites++;

### Extract characters (nucleotides) carried by the two specified sequences at this column
	my $s1_char = substr($seq_hash_local->{$sample1},$looksite,1);
	my $s2_char = substr($seq_hash_local->{$sample2},$looksite,1);

### Convert both characters to uppercase
	$s1_char = uc($s1_char);
	$s2_char = uc($s2_char);

	my $chars_column = $s1_char.$s2_char;

	if ($chars_column =~ /[^atgc]/i)
	{
	    die "\nError in the subroutine get_pairwise_distance: one or both provided sequences (IDs $sample1 and $sample2) possess non-canonical symbols at the site: $looksite.\n Stopped.\n";
	} 

	unless ($chars_column =~ /^[atgc]{2}$/i)
	{
	    die "\nError in the subroutine get_pairwise_distance: site with the index $looksite does not conform to the expected format. Stopped.\n";
	} 

### If the two analyzed sequences carry different nucleotides at the currently processed column, increment the counter for the total number of nucleotide differences
	if ($s1_char ne $s2_char)
	{
	    # print "$looksite $s1_char $s2_char\n";	    
	    $ndiverged_sites++;
	}

    }

### Compute proportion of nucleotide differences between the two sequences
    my $diverged_fraction = 0;

    eval { $diverged_fraction = $ndiverged_sites / $nlookedup_sites; }; warn $@ if $@;

    # print "nlookedup_sites $nlookedup_sites\n";
    # print "ndiverged_sites $ndiverged_sites\n";
    # print "diverged_fraction $diverged_fraction\n";

    return ($diverged_fraction,$ndiverged_sites,$nlookedup_sites);

}


