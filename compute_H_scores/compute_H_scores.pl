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
use Cwd;
use Data::Dumper;

my ($new_vcf) = @ARGV;

unless (defined($new_vcf))
{
    die "\nInput file is not set. Stopped.";
}

print "\n";
print "#######################################################################################################\n";
print "Run compute_H_scores.pl\n";
print "\n";
print "#######################################################################################################\n";
print "#######################################################################################################\n";

print "Compute H-scores using triallelic sites carrying all three heterozygous genotypes ('0/1','0/2','1/2')\namong analyzed individuals.\n";

print "\nThis script searches an input VCF file for triallelic sites carrying all three heterozygous genotypes\nthat harbor only one private heterozygous genotype.\n";
print "That is, such sites where three heterozygotes are found among the analyzed individuals and the least\nfrequent of the three heterozygous genotypes is present in a single individual with the next frequent\ngenotype present at least in two individuals.\n";
print "This subset of sites is used to compute H-scores.\n";


print "\n";
print "For further information, see README.md\n";
print "#######################################################################################################\n";
print "#######################################################################################################\n";
print "\n";
print "Started\n";
system("date");
print "\n";
print "#######################################################################################################\n";

print "\n";
print "#######################################################################################################\n";
print "Input VCF file to process:\n$new_vcf\n";
print "\n";

my $cur = cwd();

print "Current working directory:\n$cur\n";


#######################################################################################################
### Specify output file names

### Name of the output file to report per sample numbers of sites with the least frequent heterozygous genotype private to this sample
my $output_samples_least_frequent_heterozygote_counts = "$new_vcf.private.het.counts.txt";

# print "\n";
# print "output_samples_least_frequent_heterozygote_counts: $output_samples_least_frequent_heterozygote_counts\n";
open (HETUNI,">$output_samples_least_frequent_heterozygote_counts") or die $!;

### Name of the output file to print out H-scores (normalized by the total number of sites with three HETs where only one private heterozygous genotype exists)
my $output_pairs_contributors_combinations_weighted_counts_normalized_list="$new_vcf.H.scores.txt";

# print "\n";
# print "output_pairs_contributors_combinations_weighted_counts_normalized_list: $output_pairs_contributors_combinations_weighted_counts_normalized_list\n";
open (HETCONTRWN_LIST,">$output_pairs_contributors_combinations_weighted_counts_normalized_list") or die $!;

print HETCONTRWN_LIST "SAMPLE1\tSAMPLE2\tH-SCORE\n";

### Use $separator to extract the genotype (GT) from genotype fields
### According to the VCF specification, GT is the first sub-field of the genotype field
my $separator = ':';

### Retrieve a line from the VCF header containing sample IDs 
my $format_line = `grep '#CHROM' $new_vcf`;
chomp($format_line);

my @allowed_genotypes_extended = ('0/0', '0/1', '1/1', '0/2', '2/2', '1/2');

my @allowed_heterozygotes = ('0/1','0/2','1/2');
my @allowed_homozygotes = ('0/0','1/1','2/2');

print "\n";
print "#######################################################################################################\n";
print "Allowed types of genotypes: @allowed_genotypes_extended\n";
print "Genotypes other than these are not supposed to be present in the input VCF file\n";
print "#######################################################################################################\n";
print "\n";

my $allowed_genotypes_hash;

foreach my $allowed_genotype0 (@allowed_genotypes_extended)
{
    $allowed_genotypes_hash->{$allowed_genotype0} = 1;
}

#######################################################################################################
### For each allowed type of genotypes, calculate its weight and store it in a hash
### These data are used further to compute probabilities of individuals to have contributed alleles to a private heterozygous genotype### Weights are computed using the 'get_genotype_weight' subroutine

my $genotypes_weights = { map { 

    $_ => get_genotype_weight($_); 

} @allowed_genotypes_extended };


# print "\n\n";
# print Dumper($genotypes_weights);

#######################################################################################################
### For each possible least frequent heterozygous genotype private to a single individual, determine which genotypes (other than that) could potentially contribute alleles to this genotype.
### Only triallelic sites carrying three HETs are considered here.
### For example, let us assume that the '1/2' genotype is the only private heterozygote among the three found at a triallelic site.
### This implies that '1/2' genotype is present only in a single individual. In this case, the first allele of '1/2' ('1') could have been potentially contributed by samples carrying genotypes '1/1' or '0/1', and the second allele ('2') could have been contributed by samples carrying genotypes '2/2' or '0/2'. 
### Store the corresponding information in a hash.

my $heterozygous_genotypes_possible_contributors;

$heterozygous_genotypes_possible_contributors->{'0/1'}->{'0'} = ['0/0', '0/2'];
$heterozygous_genotypes_possible_contributors->{'0/1'}->{'1'} = ['1/1', '1/2'];

$heterozygous_genotypes_possible_contributors->{'0/2'}->{'0'} = ['0/0', '0/1'];
$heterozygous_genotypes_possible_contributors->{'0/2'}->{'2'} = ['2/2', '1/2'];

$heterozygous_genotypes_possible_contributors->{'1/2'}->{'1'} = ['1/1', '0/1'];
$heterozygous_genotypes_possible_contributors->{'1/2'}->{'2'} = ['2/2', '0/2'];

# print "\n\n";
# print Dumper($heterozygous_genotypes_possible_contributors);

### Retrieve a list of samples present in the input VCF file
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@samples_new) = split/\s+/,$format_line;

### Number of samples present in the input VCF file
my $dim_samples = @samples_new;

print "IDs of samples contained in the input VCF file:\n@samples_new\n";
print "\nNumber of samples in the input VCF file: $dim_samples\n\n";


### Map indices of samples in the VCF file to sample IDs
my $samples_to_ind_map;
my $ind_to_sample_map;

for (my $sind=0;$sind<@samples_new;$sind++)
{
    my $ind_sample = $samples_new[$sind];
    $samples_to_ind_map->{$ind_sample} = $sind;
    $ind_to_sample_map->{$sind} = $ind_sample;
}


#######################################################################################################
### Read the input VCF file with triallelic sites carrying all three heterozygous genotypes

### Total number of VCF entries
my $genotype_lines = 0;

### Total number of VCF entries with >2 alternate alleles (terallelic sites) - these will be skipped, if present
my $genotype_lines_to_skip = 0; 

### Total number of sites carrying three heterozygotes
my $ngenotypes_three_heterozygotes = 0;

### Total number of sites carrying three heterozygotes, where exists only one private heterozygous genotype 
### That is, such sites with three heterozygotes, where the least frequent of the three heterozygotes is present in a single individual with the next frequent heterozygous genotype present at least in two individuals
my $ngenotypes_three_heterozygotes_single_unique = 0; 


### Reference to a hash storing per sample numbers of sites with 3 HETs, such that the least frequent heterozygous genotype is private to this sample
my $samples_three_heterozygotes_leastfrequent_counts; 

### Reference to a hash storing raw pairwise H-scores (not normalized by the number of sites with 3 HETs)
my $hetal12_pairs_counts_weighted;


open (VCF,"$new_vcf") or die $!;

GENO: while (<VCF>)
{
    chomp();

    unless (/^#/)
    {
	$genotype_lines++;

	my ($CHROM1,$POS1,$ID1,$REF1,$ALT1,$QUAL1,$FILTER1,$INFO1,$FORMAT1,@genotypes) = split/\s+/;

	my $dim_genotypes = scalar(@genotypes);

	unless ($dim_genotypes == $dim_samples)
	{
	    die "\nError: number of genotype calls at genotype line $genotype_lines ($CHROM1 $POS1) isn't equal to the number of samples: $dim_genotypes vs $dim_samples. Stopped.";
	}

	my @ALTVARS = split/\,/, $ALT1;

### Number of alternative variants found for the current site - is expected to be equal to 2 for triallelic sites
	my $ALT_DIM = scalar(@ALTVARS);

### Skip tetraallelic sites
	unless ($ALT_DIM <= 2)
	{
	    print "\nmultiallelic site ($CHROM1 $POS1) with more than 2 ALT variants: $REF1 -> $ALT1, skip..\n";
	    $genotype_lines_to_skip++;
	    next GENO;
	}

### Reference to a hash storing information on genotypes found at this site
	my $genotypes_hash_local;

### Reference to a hash storing information on genotype counts for this site
	my $genotypes_counts_local;

### Check that all genotypes at $CHROM1 $POS1 are among allowed ones	
### If any of the genotypes is not found among allowed ones, skip the site
### Remember which genotypes are present at this site
      GINDS: for (my $gind=0;$gind<@genotypes;$gind++)
      { 
	  my $new_geno = $genotypes[$gind];

	  my $new_sample_for_ind = $ind_to_sample_map->{$gind};

	  my @format_fields = split/$separator/,$new_geno;
	  my $new_geno_short = $format_fields[0];

	  unless (exists($allowed_genotypes_hash->{$new_geno_short}))
	  {
	      next GENO;
	  }

	  $genotypes_hash_local->{$new_geno_short}->{$new_sample_for_ind} = 1;
	  $genotypes_counts_local->{$new_geno_short}++;  
      }

### Process a site further only if it carries all three heterozygotes (encoded as '0/1', '0/2' and '1/2')
	if ((exists($genotypes_hash_local->{'0/1'})) && (exists($genotypes_hash_local->{'0/2'})) && (exists($genotypes_hash_local->{'1/2'})))
	{
	    $ngenotypes_three_heterozygotes++;

### Retrieve frequencies of three heterozygous genotypes at this site 
	    my @heterozygotes_counts_local = map { $genotypes_counts_local->{$_}; } @allowed_heterozygotes;

### Sort the three heterozygous genotypes according to their frequencies at this site
	    my @frequency_sorted_heterozygotes = sort { ($genotypes_counts_local->{$a} <=> $genotypes_counts_local->{$b}) || ($a cmp $b) } @allowed_heterozygotes;

### Determine which heterozygote is least frequent among the three 
	    my $least_frequent_heterozygote = $frequency_sorted_heterozygotes[0];

### Genotype count for the least frequent heterozygote
	    my $least_frequent_heterozygote_count = $genotypes_counts_local->{$least_frequent_heterozygote};

### The next frequent heterozygote
	    my $next_least_frequent_heterozygote = $frequency_sorted_heterozygotes[1];

### Genotype count for the next frequent heterozygote
	    my $next_least_frequent_heterozygote_count = $genotypes_counts_local->{$next_least_frequent_heterozygote};

	    my $most_frequent_heterozygote = $frequency_sorted_heterozygotes[2];
	    my $most_frequent_heterozygote_count = $genotypes_counts_local->{$most_frequent_heterozygote};

	    # print "#############################################################################################################\n";
	    # print "$_\n";
	    # print "@allowed_heterozygotes -> @heterozygotes_counts_local\n";
	    # print "frequency_sorted_heterozygotes: @frequency_sorted_heterozygotes\n";
	    # print "least_frequent_heterozygote $least_frequent_heterozygote; least_frequent_heterozygote_count $least_frequent_heterozygote_count\n";
	    # print "next_least_frequent_heterozygote $next_least_frequent_heterozygote; next_least_frequent_heterozygote_count $next_least_frequent_heterozygote_count\n";
	    # print "\n";

### Process a site further only	if the least frequent of the three heterozygous genotypes is present in a single individual with the next frequent genotype present at least in two individuals    	   	    
	    if (($least_frequent_heterozygote_count == 1) && ($least_frequent_heterozygote_count < $next_least_frequent_heterozygote_count))
	    {
		$ngenotypes_three_heterozygotes_single_unique++;

		my @samples_least_frequent_heterozygote = keys(%{$genotypes_hash_local->{$least_frequent_heterozygote}});

### Recheck that the least frequent of the three heterozygous genotypes is indeed found only in a single individual	     
		unless (scalar(@samples_least_frequent_heterozygote) == 1)
		{
		    die "\nSize of the array containing IDs of samples presumably carrying a private heterozygous genotype ($least_frequent_heterozygote) is not equal to 1 at $CHROM1 $POS1:\n@samples_least_frequent_heterozygote.\nStopped.\n";
		}

### Increment the counter of sites with the least frequent heterozygous genotype private to the individual carrying such genotype at this site
		foreach (@samples_least_frequent_heterozygote)
		{
		    $samples_three_heterozygotes_leastfrequent_counts->{$_}++;
		}

		my $het_allele1;
		my $het_allele2;

### Check that the the least frequent heterozygous genotype conforms to the expected format
### Extract the two alleles constituting the least frequent heterozygous genotype at this site 		
		if ($least_frequent_heterozygote =~ /^([012])\/([012])$/)
		{
		    $het_allele1 = $1;
		    $het_allele2 = $2;

		    if ($het_allele1 eq $het_allele2)
		    {
			die "\nIndividual components of a presumably heterozygous genotype ($least_frequent_heterozygote) at $CHROM1 $POS1 are identical: $het_allele1 vs $het_allele2. Stopped.";
		    }

		}
		else
		{
		    die "\nPresumably heterozygous genotype ($least_frequent_heterozygote) at $CHROM1 $POS1 does not conform to the expected format. Stopped.";
		}

### Declare arrays containing genotypes that in principle could contribute the first and the second allele to the private heterozygote respectively 
### Some of these genotypes may in fact be missing from the pool of genotypes for the current site
### The data are retrieved from the hash populated at the beginning of the script (hash reference $heterozygous_genotypes_possible_contributors)
		my @allele1_genotypes_possible_contributors;
		my @allele2_genotypes_possible_contributors;

		if ((exists($heterozygous_genotypes_possible_contributors->{$least_frequent_heterozygote}->{$het_allele1})) && (exists($heterozygous_genotypes_possible_contributors->{$least_frequent_heterozygote}->{$het_allele2})))
		{
		    @allele1_genotypes_possible_contributors = @{$heterozygous_genotypes_possible_contributors->{$least_frequent_heterozygote}->{$het_allele1}};		
		    @allele2_genotypes_possible_contributors = @{$heterozygous_genotypes_possible_contributors->{$least_frequent_heterozygote}->{$het_allele2}};		
		}
		else
		{
		    die "\nThere are no predefined genotypes that could contribute alleles to the private heterozygous genotype $least_frequent_heterozygote ($CHROM1:$POS1).\nCheck the \$heterozygous_genotypes_possible_contributors hash reference. Stopped.";
		}


### Array used to store IDs of samples carrying any of genotypes that could have contributed the first allele to the private heterozygous genotype found at this site
		my @combined_allele1_possible_contributor_samples;

### Reference to a hash used to store weights for individual samples carrying genotypes that could have contributed the first allele to the private heterozygous genotype found at this site
### All samples with the same genotype obtain the same weight
### Samples with genotypes absent from the list of potential 'contributors' are not considered (effectively assigned 0 weight)
		my $allele1_possible_contributor_samples_weights;

### Variable used to store the overall weight for all samples that could have contributed the first allele to the private heterozygous genotype found at this site
### The value of this variable is used further to normalize weights for individual samples 
### Overall weight is computed by summing up weights assigned to individual samples that could have contributed the first allele to the private heterozygous genotype found at this site
		my $totalsum_allele1_possible_contributor_samples_weights = 0;

### Iterate over genotypes that could have contributed the first allele of the private heterozygous genotype found at this site 
		foreach my $new_allele1_possible_contributor (@allele1_genotypes_possible_contributors)
		{

### Retrieve the weight assigned to this genotype (1 for HETs, 2 for HOMs)
### The data are retrieved from the hash populated at the beginning of the script (hash reference $genotypes_weights)
		    my $new_allele1_possible_contributor_weight = $genotypes_weights->{$new_allele1_possible_contributor};

### Retrieve IDs of samples carrying this genotype at the current site, if any 
		    my @new_allele1_possible_contributor_samples = keys(%{$genotypes_hash_local->{$new_allele1_possible_contributor}});

### Iterate over samples carrying this genotype at the current site
### Remember the weight assigned to the current sample and increment the overall weight 
### The overall weight is used further to normalize weigths assigned to individual samples 
### All samples with the same genotype get the same weight

		    if (scalar(@new_allele1_possible_contributor_samples > 0))
		    {
			push @combined_allele1_possible_contributor_samples, @new_allele1_possible_contributor_samples;

			foreach (@new_allele1_possible_contributor_samples)
			{
			    $allele1_possible_contributor_samples_weights->{$_} = $new_allele1_possible_contributor_weight;
			    $totalsum_allele1_possible_contributor_samples_weights += $new_allele1_possible_contributor_weight;
			}

		    }

		}


### Array used to store IDs of samples carrying any of genotypes that could have contributed the second allele to the private heterozygous genotype found at this site
		my @combined_allele2_possible_contributor_samples;

### Reference to a hash used to store weights for individual samples carrying genotypes that could have contributed the second allele to the private heterozygous genotype found at this site
### All samples with the same genotype obtain the same weight
### Samples with genotypes absent from the list of potential 'contributors' are not considered (effectively assigned 0 weight)
		my $allele2_possible_contributor_samples_weights;

### Variable used to store the overall weight for all samples that could have contributed the second allele to the private heterozygous genotype found at this site
### The value of this variable is used further to normalize weights for individual samples 
### Overall weight is computed by summing up weights assigned to individual samples that could have contributed the second allele to the private heterozygous genotype found at this site
		my $totalsum_allele2_possible_contributor_samples_weights = 0;

### Iterate over genotypes that could have contributed the second allele of the private heterozygous genotype found at this site
		foreach my $new_allele2_possible_contributor (@allele2_genotypes_possible_contributors)
		{

### Retrieve the weight assigned to this genotype (1 for HETs, 2 for HOMs)
### The data are retrieved from the hash populated at the beginning of the script (hash reference $genotypes_weights)
		    my $new_allele2_possible_contributor_weight = $genotypes_weights->{$new_allele2_possible_contributor};

### Retrieve IDs of samples carrying this genotype at the current site, if any
		    my @new_allele2_possible_contributor_samples = keys(%{$genotypes_hash_local->{$new_allele2_possible_contributor}});


### Iterate over samples carrying this genotype at the current site
### Remember the weight assigned to the current sample and increment the overall weight
### The overall weight is used further to normalize weigths assigned to individual samples
### All samples with the same genotype get the same weight

		    if (scalar(@new_allele2_possible_contributor_samples>0))
		    {
			push @combined_allele2_possible_contributor_samples, @new_allele2_possible_contributor_samples;

			foreach (@new_allele2_possible_contributor_samples)
			{
			    $allele2_possible_contributor_samples_weights->{$_} = $new_allele2_possible_contributor_weight;
			    $totalsum_allele2_possible_contributor_samples_weights += $new_allele2_possible_contributor_weight;
			}

		    }

		}

		### Number of samples that could have possibly contributed allele1 to the given heterozygous genotype
		my $dim_samples_het_allele1 = scalar(@combined_allele1_possible_contributor_samples);
 
		### Number of samples that could have possibly contributed allele2 to the given heterozygous genotype
		my $dim_samples_het_allele2 = scalar(@combined_allele2_possible_contributor_samples);
		

		unless (($dim_samples_het_allele1 > 0) && ($dim_samples_het_allele2 > 0))
		{

		    die "\nThere are no samples that could have possibly contributed to one or both components of the heterozygote: het_allele1 $het_allele1 -> combined_allele1_possible_contributor_samples (@combined_allele1_possible_contributor_samples), het_allele2 $het_allele2 -> combined_allele2_possible_contributor_samples (@combined_allele2_possible_contributor_samples), die..";
		}

		
		my $total_al12_pairs_combinations = $dim_samples_het_allele1*$dim_samples_het_allele2;

		my $total_al12_pairs_comparisons = 0;
		
		my $totalsum_pairs_weights = 0;

### Iterate over all possible pairs of samples SAMPLE1-SAMPLE2, where:
### - SAMPLE1 is one of the samples that could have potentially contributed the first allele to the private heterozygous genotype found at this site
### - SAMPLE2 is one of the samples that could have potentially contributed the second allele to the private heterozygous genotype found at this site
### Pairs of samples other than these are ignored, effectively getting score 0

		foreach my $hetal1_sample (@combined_allele1_possible_contributor_samples)
		{

### Compute the probability for SAMPLE1 to have contributed the first allele to the private heterozygous genotype 
### The probability is computed as the weight assigned to SAMPLE1 normalized the sum of weights for all samples that could contribute the first allele
		    my $hetal1_sample_weight = $allele1_possible_contributor_samples_weights->{$hetal1_sample};

		    my $hetal1_sample_weight_normalized;

		    eval { $hetal1_sample_weight_normalized = $hetal1_sample_weight / $totalsum_allele1_possible_contributor_samples_weights; }; warn $@ if $@;

		    
		    foreach my $hetal2_sample (@combined_allele2_possible_contributor_samples)
		    {

### SAMPLE1 should be different from SAMPLE2, as the only genotype that can carry both the first and the second allele constituting the private heterozygous genotype is the one corresponding to the private heterozygote itself
			if ($hetal1_sample eq $hetal2_sample)
			{

			    die "\nSample $hetal1_sample is listed among the samples that could have contributed the first allele to the private heterozygote $least_frequent_heterozygote as well as among the samples that could have contributed the second allele at $CHROM1 $POS1:\nSamples that could contribute the first allele ($het_allele1): @combined_allele1_possible_contributor_samples\nSamples that could contribute the second allele ($het_allele2): @combined_allele2_possible_contributor_samples.\nStopped.\n";
			}

			$total_al12_pairs_comparisons++;

### Compute the probability for SAMPLE2 to have contributed the second allele to the private heterozygous genotype 
### The probability is computed as the weight assigned to SAMPLE2 normalized the sum of weights for all samples that could contribute the second allele
			my $hetal2_sample_weight = $allele2_possible_contributor_samples_weights->{$hetal2_sample};

			my $hetal2_sample_weight_normalized;

			eval { $hetal2_sample_weight_normalized = $hetal2_sample_weight / $totalsum_allele2_possible_contributor_samples_weights; }; warn $@ if $@;


			$totalsum_pairs_weights+=($hetal1_sample_weight_normalized*$hetal2_sample_weight_normalized);

### Compute the probability that the private heterozygote at this site is a result of contamination/exchange involving SAMPLE1 and SAMPLE2
			my $al12_pair_weighted_score = $hetal1_sample_weight_normalized*$hetal2_sample_weight_normalized;
			
### Increment the overall score for SAMPLE1-SAMPLE2 (and symmetrically, SAMPLE2-SAMPLE1) by the computed probability
			$hetal12_pairs_counts_weighted->{$hetal1_sample}->{$hetal2_sample}+=$al12_pair_weighted_score;
			$hetal12_pairs_counts_weighted->{$hetal2_sample}->{$hetal1_sample}+=$al12_pair_weighted_score;
								    
		    }


		}

	    }

	}
	
	else
	{
	    next GENO;
	}

    }

}

close(VCF);

print "#######################################################################################################\n";
print "\n";
print "Total number of VCF entries (data lines with genotype information): $genotype_lines\n";
print "Total number of VCF entries corresponding to sites with >2 ALT alleles (these were skipped, if present): $genotype_lines_to_skip\n";

print "\n";
print "Total number of sites carrying three heterozygotes: $ngenotypes_three_heterozygotes\n";
print "Total number of sites carrying three heterozygotes where only one private heterozygous genotype exists (these were actually processed): $ngenotypes_three_heterozygotes_single_unique\n";
print "\n";

print "#######################################################################################################\n";
print "Per sample numbers of sites with the least frequent heterozygous genotype private to this sample (out of $ngenotypes_three_heterozygotes_single_unique):\n";
print Dumper($samples_three_heterozygotes_leastfrequent_counts);

print "#######################################################################################################\n";
### Print out per sample numbers of sites with the least frequent heterozygous genotype private to this sample to the output file
foreach my $sample_toprint (@samples_new)
{
    if (exists($samples_three_heterozygotes_leastfrequent_counts->{$sample_toprint}))
    {
	print HETUNI "$sample_toprint\t$samples_three_heterozygotes_leastfrequent_counts->{$sample_toprint}\n";
    }   
    else
    {
	print HETUNI "$sample_toprint\t0\n";
    }

}

print "\n";
print "#######################################################################################################\n";
print "Per sample numbers of sites with three heterozygotes and the least frequent heterozygous genotype private to this sample were written to the file:\n";
print "$output_samples_least_frequent_heterozygote_counts\n";


if ($ngenotypes_three_heterozygotes_single_unique < 50)
{
    die "\nToo few triallelic sites carrying all three heterozygous genotypes and harboring only one private heterozygous genotype ($ngenotypes_three_heterozygotes_single_unique) to compute H-scores.\nStopped.\n";
}

#######################################################################################################
### Compute normalized H-scores for each pair of samples, dividing raw H-scores by the total number of sites with three HETs where only one private least frequent heterozygous genotype exists

my $total_sum_unique_pairwise_weights = 0;
my $total_sum_unique_pairwise_weights_normalized = 0;

for (my $sindnew=0;$sindnew<@samples_new-1;$sindnew++)
{
    my $pair_sample1 = $samples_new[$sindnew];
    
    for (my $sindnew1=$sindnew+1;$sindnew1<@samples_new;$sindnew1++)
    {
	my $pair_sample2 = $samples_new[$sindnew1];

	my $pair_al12_total_contribution_weighted_count = $hetal12_pairs_counts_weighted->{$pair_sample1}->{$pair_sample2};

### If H-score is not set, set H-score to 0 
	unless ($pair_al12_total_contribution_weighted_count)
	{
	    $pair_al12_total_contribution_weighted_count = 0;
	}

	my $pair_al12_total_contribution_weighted_count_normalized; 

	eval { $pair_al12_total_contribution_weighted_count_normalized = $pair_al12_total_contribution_weighted_count / $ngenotypes_three_heterozygotes_single_unique; }; warn $@ if $@;

	$total_sum_unique_pairwise_weights += $pair_al12_total_contribution_weighted_count;
	$total_sum_unique_pairwise_weights_normalized += $pair_al12_total_contribution_weighted_count_normalized;

	print HETCONTRWN_LIST "$pair_sample1\t$pair_sample2\t$pair_al12_total_contribution_weighted_count_normalized\n"

    }    

}

close(HETUNI);
close(HETCONTRWN_LIST);

print "\n";
print "#######################################################################################################\n";
print "Pairwise H-scores were written to the file:\n";
print "$output_pairs_contributors_combinations_weighted_counts_normalized_list\n";


# print "\n";
# print "##############################################################################\n";
# print "Check that the total sum of normalized pairwise H-scores is equal to 1:\n";
# print "total_sum_unique_pairwise_weights_normalized: $total_sum_unique_pairwise_weights_normalized\n";
# print "##############################################################################\n\n";


print "\n";
print "#######################################################################################################\n";
print "\nFinished\n";
system("date");



### Subroutines

sub get_genotype_weight
{
    my ($newgenotype_input) = @_;

    unless ($newgenotype_input =~ /^([012])\/([012])$/)
    {
	die "\nError in the 'get_genotype_weight' subroutine: genotype $newgenotype_input does not conform to the expected format. Stopped.";
    }

    my ($gsplit1,$gsplit2) = split/\//,$newgenotype_input;
 
    my $gweight;

### Assign the value 2, if a genotype is homozygous
    if ($gsplit1 eq $gsplit2)
    {
	$gweight = 2;
    }
### Assign the value 1, otherwise
    else
    {
	$gweight = 1;
    }

    return($gweight);
}
