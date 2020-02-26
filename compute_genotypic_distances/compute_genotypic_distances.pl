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
use List::Util qw(sum min max);

print "\n";
print "########################################################################################################################\n";
print "########################################################################################################################\n";
print "Compute genotypic distances for all pairs of samples present in the input VCF file\n";
print "\n";
print "This script is supposed to work only with biallelic SNP sites and monomorphic sites genotyped in all the samples contained in the VCF file\n";
print "INDELs are not supposed to be present in the input VCF file\n";
print "Only a single entry is supposed to be present for each included genomic site\n";
print "\n";

print "Sample names must not contain dots (.)\n";
print "########################################################################################################################\n";
print "########################################################################################################################\n";

#### This script is supposed to work with sites genotyped in all the samples present in the input VCF file
#### Sites with missing genotypes should be filtered out before running the script
#### Only biallelic SNP sites and monomorphic sites are expected to be present in the input VCF file
#### All VCF entries should correspond to unique genomic sites
#### Multiallelic and INDEL sites should be filtered out before running the script

my ($new_vcf) = @ARGV;

print "\n";
print "Started\n";
system("date");

print "\n";
print "########################################################################################################################\n";
print "Input VCF file to process: $new_vcf\n";
print "\n";

my $cur = cwd();

print "Current working directory: $cur\n";

my $format_line=`grep '^#CHROM' $new_vcf`;
chomp($format_line);

my @allowed_genotypes=('0/0', '0/1', '1/1');

print "\n";
print "########################################################################################################################\n";
print "Allowed types of genotypes: @allowed_genotypes\n";
print "Genotypes other than these are not supposed to be present in the input VCF file\n";
print "########################################################################################################################\n";
print "\n";

my $allowed_genotypes_hash;

foreach my $allowed_genotype0 (@allowed_genotypes)
{
    $allowed_genotypes_hash->{$allowed_genotype0} = 1;
}

my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@samples_new) = split/\s+/,$format_line;

# print "format_line: $format_line\n";
print "Names of the samples contained in the input VCF file: @samples_new\n";

my $dim_samples = @samples_new;
print "\nNumber of samples in the input VCF file: $dim_samples\n\n";

my $samples_to_ind_map;
my $ind_to_sample_map;

for (my $sind=0;$sind<@samples_new;$sind++)
{
    my $ind_sample = $samples_new[$sind];
    # print "$sind -> $ind_sample\n";
    $samples_to_ind_map->{$ind_sample} = $sind;
    $ind_to_sample_map->{$sind} = $ind_sample;
}

open (VCF,"$new_vcf") or die "Can't open $new_vcf: $!";

my $genotype_lines = 0;
my $multiallelic_or_indel_genotype_lines = 0;

my $used_genotype_lines = 0;
my $used_genotype_lines_monomorphic = 0;
my $used_genotype_lines_segregating = 0;

### Reference to a hash storing pairwise genotypic distances for all possible pairs of samples
my $samples_pairwise_distances_mutual_hash;  

my $total_monomorphic_effective_sites_number = 0;
my $total_monomorphic_sites_number = 0;

### Use $separator to extract the genotype (GT) from genotype fields 
### According to the VCF specification, GT is the first sub-field of the genotype field
my $separator=':';

GENO: while (<VCF>)
{
    chomp();

###Header lines are skipped
    unless (/^#/)
    {
	$genotype_lines++;

	my ($CHROM1,$POS1,$ID1,$REF1,$ALT1,$QUAL1,$FILTER1,$INFO1,$FORMAT1,@genotypes) = split/\s+/;

	my $ref_dim = length($REF1);
	my $alt_dim = length($ALT1);

	unless (($alt_dim == 1) && ($ref_dim == 1))
	{
	    warn "\nWarning: multiallelic or INDEL site (ref_dim==$ref_dim; alt_dim==$alt_dim) $CHROM1 $POS1: $REF1 -> $ALT1, skip..\n";
	    $multiallelic_or_indel_genotype_lines++;
	    next GENO;
	}
	
	my $dim_genotypes = scalar(@genotypes);

	unless ($dim_genotypes == $dim_samples)
	{
	    die "\nError: number of genotype calls at genotype line $genotype_lines ($CHROM1 $POS1) isn't equal to the number of samples: $dim_genotypes vs $dim_samples. Stopped.";
	}

### Check that all genotypes at $CHROM1 $POS1 are among allowed ones

      GINDS: for (my $gind=0;$gind<@genotypes;$gind++)
      {
          my $new_geno = $genotypes[$gind];

          my $new_sample_for_ind = $ind_to_sample_map->{$gind};

          my @genotype_fields = split/$separator/,$new_geno;
          my $new_geno_short = $genotype_fields[0];

          unless (exists($allowed_genotypes_hash->{$new_geno_short}))
          {
	      die "\nError: unexpected genotype at $CHROM1 $POS1: $new_sample_for_ind $new_geno_short ($new_geno). Stopped.";
          }

      }

	$used_genotype_lines++;

### Calculate pairwise distances for all possible pairs of samples among those present in the VCF file 
### Whole-genome distances are computed by summing the distances computed for individual sites 
### Monomorphic (reference monomorphic, i.e. with genotype 0/0 in all of the samples) and polymorphic sites are processed separately

	if ($ALT1 eq '.')
	{

### Process monomorphic reference sites
### Check whether a site is truly monomorphic (i.e. contains only 0/0 genotypes)

	    my $site_genotypes_index;

### Iterate over all samples and extract their genotypes
	  GINDS0: for (my $gind0=0;$gind0<@genotypes;$gind0++)
	  {
	      my $new_geno0 = $genotypes[$gind0];
	      
	      my @genotype_fields0 = split/$separator/,$new_geno0;

###Variable $new_geno_short0 stores the genotype (GT) for the processed sample
	      my $new_geno_short0 = $genotype_fields0[0];

	      $site_genotypes_index->{$new_geno_short0} = 1;

	  }

	    my @site_unique_genotypes = keys(%{$site_genotypes_index});
	    my $site_unique_genotypes_string = join('.', @site_unique_genotypes);

	    unless ($site_unique_genotypes_string eq '0/0')
	    {
		die "\nError: unexpected genotypes at a presumably monomorphic site $CHROM1 $POS1: @site_unique_genotypes, ($_). Stopped.";
	    }
	    else
	    {
		$used_genotype_lines_monomorphic++;
	    
		$total_monomorphic_sites_number++;

### Effective number of monomorphic sites is incremented by 2, as each VCF entry corresponds to a diploid genotype call
### Therefore, effectively two sites (one from each of the two haplotypes are compared)
### Only the effective number of sites is incremented here, distance is not incremented as at a monomoprhic reference site all the samples have the same genotype
		$total_monomorphic_effective_sites_number+=2;
	    }

	}
	elsif ($ALT1 =~ /^[ATGC]$/i)
	{

	    $used_genotype_lines_segregating++;

### Iterate over all different pairs of samples contained in the input VCF file
### Compute all pairwise genotypic distances at the considered site 
### The pairwise genotypic distance is computed as the difference between the numbers of non-reference variants harbored by these two samples at this site 

	  GINDS1: for (my $gind1=0;$gind1<@genotypes;$gind1++)
	  {
	      my $new_geno1 = $genotypes[$gind1];
	  
### Variable $new_sample_for_ind1 stores the name for the first sample in the pair
	      my $new_sample_for_ind1 = $ind_to_sample_map->{$gind1};

	      my @genotype_fields1 = split/$separator/,$new_geno1;

### Variable $new_geno_short1 stores the genotype (GT) for the first sample in the pair
	      my $new_geno_short1 = $genotype_fields1[0];

	    GINDS2: for (my $gind2=$gind1;$gind2<@genotypes;$gind2++)
	    {
		my $new_geno2 = $genotypes[$gind2];
	  
### Variable $new_sample_for_ind2 stores the name for the second sample in the pair
		my $new_sample_for_ind2 = $ind_to_sample_map->{$gind2};

		my @genotype_fields2 = split/$separator/,$new_geno2;

### Variable $new_geno_short2 stores the genotype (GT) for the second sample in the pair
		my $new_geno_short2 = $genotype_fields2[0];

### Distance between the two genotypes is computed using subroutine 'genotypes_compare'
		my ($pair_effective_sites_number12,$pair_gntp_distance12) = genotypes_compare($new_geno_short1, $new_geno_short2);

		unless ($pair_effective_sites_number12 == 2)
		{
		    die "\nError: subroutine 'genotypes_compare' yields an unexpected effective number of sites at $CHROM1 $POS1: $pair_effective_sites_number12. Stopped.\n";
		}

### Check that comparison of a genotype vs itself yields a zero distance
		if (($gind1 == $gind2) && ($pair_gntp_distance12 != 0))
		{
		    die "\nError: comparison of a genotype vs itself yields a non-zero distance at $CHROM1 $POS1: $new_sample_for_ind1 vs $new_sample_for_ind2; $gind1 vs $gind2; $new_geno_short1 vs $new_geno_short2 gives distance = $pair_gntp_distance12. Stopped.";
		}

### Add a distance between the two samples in the pair ($new_sample_for_ind1 $new_sample_for_ind2) at the processed site $CHROM1 $POS1 to the overall genotypic distance between these two samples
### The overall genotypic distance is incremented in both directions ($new_sample_for_ind1 -> $new_sample_for_ind2) and ($new_sample_for_ind2 -> $new_sample_for_ind1), so that the two symmetrical distances are stored in the hash

		my $new_pair_direct="$new_sample_for_ind1.$new_sample_for_ind2";
		my $new_pair_reciprocal="$new_sample_for_ind2.$new_sample_for_ind1";

		my $distinct_pairs_href = { map { $_ => 1 } ($new_pair_direct, $new_pair_reciprocal) };
		my @distinct_pairs_array = keys(%{$distinct_pairs_href});


		foreach my $new_pair_mutual (@distinct_pairs_array)
		{
		    if ($new_pair_mutual =~ /^(\S+)\.(\S+)$/)
		    {
			my ($new_sample_mutual1, $new_sample_mutual2) = ($1,$2);

### Increment the effective number of compared polymorphic sites for the given pair of samples
			$samples_pairwise_distances_mutual_hash->{$new_sample_mutual1}->{$new_sample_mutual2}->{'effective_sites_number'}+=$pair_effective_sites_number12;

### Increment the overall genotypic distance for the given pair of samples by the computed value (can be 0, 1 or 2)
			$samples_pairwise_distances_mutual_hash->{$new_sample_mutual1}->{$new_sample_mutual2}->{'genotype_distance'}+=$pair_gntp_distance12;

			$samples_pairwise_distances_mutual_hash->{$new_sample_mutual1}->{$new_sample_mutual2}->{'sites_compared'}++;

		    }
		    else
		    {
			die "\nError: ID for the pair of samples $new_pair_mutual does not conform to the expected format. Stopped.";
		    }
		}
	    }
	      
	  }
	    
	}
	else 
	{
	    die "\nError: unexpected ALT FIELD at $CHROM1 $POS1: $ALT1 ($_). Stopped.";
	}

    }

}

close(VCF);

print "########################################################################################################################\n";
print "########################################################################################################################\n";
print "\n";
print "Total number of VCF entries (data lines with genotype information): $genotype_lines\n";
print "Total number of VCF entries corresponding to INDEL or multiallelic sites (these are skipped if present): $multiallelic_or_indel_genotype_lines\n";


print "\n";
print "########################################################################################################################\n";
print "Total number of analyzed VCF entries: $used_genotype_lines\n";
print "\n";
print "Among $used_genotype_lines analyzed VCF entries:\n";
print "$used_genotype_lines_monomorphic correspond to monomorphic sites\n";
print "$used_genotype_lines_segregating correspond to polymorphic sites\n";
print "\n";


print "########################################################################################################################\n";
print "########################################################################################################################\n";
print "Report computed pairwise distances:\n";
print "\n";

### Specify file names to report pairwise genotypic distances

my $samples_string = join("\t",@samples_new);

### Report absolute genotypic distances 
my $output_pairwise_distances_file_absolute = "$new_vcf.pair.gntp.abs.dist.txt";
open (PADFA,">$output_pairwise_distances_file_absolute") or die $!;
print PADFA "\t$dim_samples\n";
print PADFA "\t$samples_string\n";

### Report genotypic distances normalized by the total effective number of analyzed sites (both monomorphic and polymorphic)
my $output_pairwise_distances_file_fractions_allsites = "$new_vcf.pair.gntp.frac.allsites.dist.txt";
open (PADFF_ALL,">$output_pairwise_distances_file_fractions_allsites") or die $!;
print PADFF_ALL "\t$dim_samples\n";
print PADFF_ALL "\t$samples_string\n";

### Report genotypic distances normalized by the total effective number of polymorphic sites
my $output_pairwise_distances_file_fractions = "$new_vcf.pair.gntp.frac.segsites.dist.txt";
open (PADFF,">$output_pairwise_distances_file_fractions") or die $!;
print PADFF "\t$dim_samples\n";
print PADFF "\t$samples_string\n";

### Index effective numbers of polymorphic sites analyzed for different pairs of samples - should be the same number for all pairs
my $pairwise_effective_sites_number_index;

### Index total effective numbers of sites (both monomorphic and polymorphic) analyzed for different pairs of samples - should be the same number for all pairs
my $pairwise_effective_sites_number_all_index;


SINDS1: for (my $sind1=0;$sind1<@samples_new;$sind1++)
{
    my $sind_sample1 = $samples_new[$sind1];

### String to store absolute genotypic distances
    my $string_to_print_distabsolute="$sind_sample1\t"; 

###  String to store distances normalized by the effective number of polymorphic sites
    my $string_to_print_distfractions="$sind_sample1\t"; 

### String to store distances normalized by the total effective number of analyzed sites, including monomorphic ones
    my $string_to_print_distfractions_all="$sind_sample1\t"; 

  SINDS2: for (my $sind2=0;$sind2<@samples_new;$sind2++)
  {
      my $sind_sample2 = $samples_new[$sind2];

### Compute overall genotypic distance normalized by the total effective number of polymorphic sites
      my $total_pairwise_genotype_distance = $samples_pairwise_distances_mutual_hash->{$sind_sample1}->{$sind_sample2}->{'genotype_distance'};
      my $total_pairwise_effective_sites_number = $samples_pairwise_distances_mutual_hash->{$sind_sample1}->{$sind_sample2}->{'effective_sites_number'};

      my $total_pairwise_genotype_distance_fraction;

      unless (defined($total_pairwise_genotype_distance))
      {
	   warn "\nWarning: absolute pairwise genotypic distance for the comparison of $sind_sample1 vs $sind_sample2 is not set!\n\n";
	  $total_pairwise_genotype_distance = 0;
      }

      if ($total_pairwise_effective_sites_number > 0)
      {
	  eval { $total_pairwise_genotype_distance_fraction = $total_pairwise_genotype_distance / $total_pairwise_effective_sites_number; }; warn $@ if $@;

      }
      else 
      {
	  $total_pairwise_effective_sites_number = 0;
	  $total_pairwise_genotype_distance_fraction='NA';
      }

### Compute genotypic distance normalized by the total effective number of sites present in the input VCF file, including monomorphic ones
      my $total_pairwise_effective_sites_number_all = $total_pairwise_effective_sites_number+$total_monomorphic_effective_sites_number;

      my $total_pairwise_genotype_distance_fraction_allsites;

      eval { $total_pairwise_genotype_distance_fraction_allsites = $total_pairwise_genotype_distance / $total_pairwise_effective_sites_number_all; }; warn $@ if $@;


      $pairwise_effective_sites_number_index->{$total_pairwise_effective_sites_number} = 1;
      $pairwise_effective_sites_number_all_index->{$total_pairwise_effective_sites_number_all} = 1;


      $string_to_print_distabsolute.="$total_pairwise_genotype_distance\t";

      $string_to_print_distfractions.="$total_pairwise_genotype_distance_fraction\t";

      $string_to_print_distfractions_all.="$total_pairwise_genotype_distance_fraction_allsites\t";

  }

    print PADFA "$string_to_print_distabsolute\n";
    print PADFF "$string_to_print_distfractions\n";
    print PADFF_ALL "$string_to_print_distfractions_all\n";

}

close(PADFA);
close(PADFF);
close(PADFF_ALL);


my @unique_pairwise_effective_sites_numbers = keys(%{$pairwise_effective_sites_number_index}); 
my @unique_pairwise_effective_sites_numbers_all = keys(%{$pairwise_effective_sites_number_all_index});

my $pairwise_effective_sites_used = 'NA';
my $pairwise_effective_sites_all_used = 'NA';

unless (scalar(@unique_pairwise_effective_sites_numbers) == 1)
{
    warn "\nWarning: different effective numbers of polymorphic sites are analyzed for different pairs of samples: @unique_pairwise_effective_sites_numbers, this should not happen..";
}
else
{
    $pairwise_effective_sites_used = $unique_pairwise_effective_sites_numbers[0];
}

### Check whether the same effective number of sites is used in all pairwise comparisons 
unless (scalar(@unique_pairwise_effective_sites_numbers_all) == 1)
{
    warn "\nWarning: different effective numbers of sites are analyzed for different pairs of samples: @unique_pairwise_effective_sites_numbers_all, this should not happen..";
}
else
{
    $pairwise_effective_sites_all_used = $unique_pairwise_effective_sites_numbers_all[0];
}

print "########################################################################################################################\n";
print "Absolute genotypic distances were written to the file:\n$output_pairwise_distances_file_absolute\n";
print "\n";

print "########################################################################################################################\n";
print "Genotypic distances normalized by the total effective number of analyzed sites (both monomorphic and polymorphic) were written to the file:\n$output_pairwise_distances_file_fractions_allsites\n";
print "\n";
print "Effective number of analyzed sites used to normalize genotypic distances: $pairwise_effective_sites_all_used\n";
print "Number of analyzed VCF entries: $used_genotype_lines\n";
print "\n";

print "########################################################################################################################\n";
print "Genotypic distances normalized by the effective number of polymorphic sites were written to the file:\n$output_pairwise_distances_file_fractions\n";
print "\n";
print "Effective number of polymorphic sites used to normalize genotypic distances: $pairwise_effective_sites_used\n";
print "Number of analyzed VCF entries corresponding to polymorphic sites: $used_genotype_lines_segregating\n";
print "\n";

print "########################################################################################################################\n";
print "\nFinished\n";
system("date");


### Subroutines

sub genotypes_compare 
{

    my ($gntp1,$gntp2)=@_;

    my $effective_sites_number;

### Computes difference in the number of non-reference calls between two genotypes
### Only diploid genotypes corresponding to biallelic SNP sites or monomorphic sites are allowed (genotypes: 0/0, 0/1, 1/1)
### VCF file is not supposed to contain phased data

    my $gntp_regex=qr/^[01]{1}\/[01]{1}$/;

    if (($gntp1 =~ /$gntp_regex/) && ($gntp2 =~ /$gntp_regex/))
    {
	$effective_sites_number = 2;
    }
    else
    {
	die "\nError in the subroutine 'genotypes_compare': genotypes do not conform to the expected format: $gntp1 vs $gntp2. Stopped.";
    }

    my @genotype_sums = map { 
	my @gntps_sepr = split/\//;
	my $gntps_sum = sum(@gntps_sepr);
	$gntps_sum;
    } ($gntp1,$gntp2);

    unless (scalar(@genotype_sums) == 2)
    {
	die "\nError in the subroutine 'genotypes_compare': unexpected number of elements in the array genotype_sums: @genotype_sums. Stopped.";
    }

    my $min_gntp_sum = min(@genotype_sums);
    my $max_gntp_sum = max(@genotype_sums);

    my $gntp_distance = $max_gntp_sum-$min_gntp_sum;

    if (($gntp_distance < 0) or ($gntp_distance > 2))
    {
	die "\nError in the subroutine 'genotypes_compare': unexpected distance between genotypes $gntp1 and $gntp2: $gntp_distance. Stopped.";
    }

    return($effective_sites_number, $gntp_distance);

}


sub mean {
    return sum(@_)/@_;
}
