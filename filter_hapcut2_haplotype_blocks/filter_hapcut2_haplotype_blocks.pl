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
use List::Util qw[min max];
use YAML;

my ($haplotype_file,$indices_conflicts_file) = @ARGV;

unless (defined($haplotype_file) && defined($indices_conflicts_file))
{
    die "\nInput files are not set. Stopped.";
}

print "\n";
print "#######################################################################################################\n";
print "Run filter_hapcut2_haplotype_blocks.pl\n";
print "#######################################################################################################\n";

print "\n";
print "This script is meant to be used for filtering of phased haplotype blocks assembled with HapCUT2.\n";

print "\nThe script takes as input two files:\n";
print "1) a haplotype file produced with HapCUT2 (https://github.com/vibansal/HapCUT2).\n";
print "2) a file containing indices of conflicting variant pairs produced with the ‘get_conflicting_variants_indices.pl’ script.\n";

print "\n";

print "The script produces two output files containing two sets of filtered haplotype blocks:\n";

print "1) output file with the extension ‘.conflicts.filtered.txt’ contains filtered Set 1.\nThis set includes those original intact haplotype blocks that were free of variants involved in conflicting variant pairs.\n";
print "\n";
print "2) output file with the extension ‘.conflicts.qual.filtered.txt’ contains filtered Set 2.\n";
print "This set includes:\na) original intact haplotype blocks that were free both of variants involved in conflicts and problematic* variants.\nb) chunks of original haplotype blocks resulting from the split of an original block at the position of a problematic variant.\n";

print "\n";
print "*A phased variant is considered problematic, if its switch or mismatch quality computed with HapCUT2 < 100. This threshold can be modified by changing the value assigned to the variable \$EXP_QUAL_VALUE inside the script.\n";

print "\n";
print "For further information, see README.md\n";

### Specify output file names to print out filtered haplotype blocks

### Name of the output file containing phased blocks filtered only based on the supplied list of conflicting variant indices
my $haplotype_file_conflicts_filtered = "$haplotype_file.conflicts.filtered.txt";
# print "\nhaplotype_file_conflicts_filtered: $haplotype_file_conflicts_filtered\n";
# print "\n";

if ($haplotype_file_conflicts_filtered eq $haplotype_file)
{
    die "\nNames of the input haplotype file and the output haplotype file filtered only for conflicts are identical: $haplotype_file_conflicts_filtered vs $haplotype_file, die..";
}

### Name of the output file containing phased blocks filtered both based on the supplied list of conflicting variant indices and based on probabilities of phasing errors computed with HapCUT2
my $haplotype_file_filtered = "$haplotype_file.conflicts.qual.filtered.txt";
# print "\nhaplotype_file_filtered: $haplotype_file_filtered\n";
# print "\n";

if ($haplotype_file_filtered eq $haplotype_file)
{
    die "\nNames of the input haplotype file and the output haplotype file filtered for conflicts and mismatch/switch qualities are identical: $haplotype_file_filtered vs $haplotype_file, die..";
}

print "#######################################################################################################\n";
print "\n";
print "Started\n";
system("date");
print "\n";
print "#######################################################################################################\n";

print "\n";
print "#######################################################################################################\n";
print "#######################################################################################################\n";
print "Read input files:\n";
print "\n";
print "#######################################################################################################\n";
print "Input file with indices of conflicting pairs of variants:\n$indices_conflicts_file\n";
print "\n";
print "Number of lines in this file (each line corresponds to a pair of conflicting variants):\n";
system("wc -l $indices_conflicts_file");

### Proceed to reading the input file with indices of conflicting pairs of variants

### Declare variables used in this section of the script

### Variable $CONFLICTING_PAIRS_INDICES_HASH is a reference to a hash storing a list of indices of conflicting variant pairs
### It is autovivified to be a hash reference, when the first element is added
### Index of a pair of conflicting variants is supposed to have been constructed by concatenating the two corresponding variant indices with a dot ('.') as a separator 
### For example, record 186.194 would correspond to a conflicting pair of variants with indices 186 and 194.
my $CONFLICTING_PAIRS_INDICES_HASH;

### Variable $CONFLICTING_SINGLE_INDICES_HASH is a reference to a hash storing a list of indices for individual variants included in any conflicting pair 
### It is autovivified to be a hash reference, when the first element is added
### For example, if the input file contained 3 records with the following indices of conflicting variant pairs: 186.194, 194.250, 270.305, the list of indices for variants included in any conflicting pair would be: 186, 194, 250, 270, 305
### Index of a variant is just its index number in the VCF file supplied to extractHAIRS to produce a fragment file further used by HapCUT2 to perform phasing
### Index of a variant is not its chromosomal position!!!
my $CONFLICTING_SINGLE_INDICES_HASH;

## Variable $CONFLICTING_SINGLE_INDICES_HASH is a reference to a hash of arrays mapping indices of individual variants to indices of conflicting variant pairs
### It is autovivified, when the first element is added
my $CONFLICTING_SINGLE_INDICES_TO_PAIRS;

### Variable $TOTAL_CONFLICTING_PAIRS_READ stores the total number of conflicting variant pairs 
my $TOTAL_CONFLICTING_PAIRS_READ = 0;

my $EXP_DIM_CONF_FIELDS = 5;

### Read input file with the list of conflicting pairs of variants
open (CONFIND,"$indices_conflicts_file") or die $!;

  while (<CONFIND>)
{
    chomp();

### Increment the total number of conflicting variant pairs
    $TOTAL_CONFLICTING_PAIRS_READ++;

    my @CONF_FIELDS = split/\;\s+/;

    my $DIM_CONF_FIELDS = scalar(@CONF_FIELDS);

### Check that the number of fields in the current line matches the expected one
    unless ($DIM_CONF_FIELDS == $EXP_DIM_CONF_FIELDS)
    {
	die "\nUnexpected number of fields in the file with the list of indices of conflicting variant pairs: $DIM_CONF_FIELDS vs $EXP_DIM_CONF_FIELDS\nat line $TOTAL_CONFLICTING_PAIRS_READ: $_. Stopped.";
    }

### Variable $CONF_IND_PAIR stores an index of the pair of conflicting variants
    my $CONF_IND_PAIR = $CONF_FIELDS[0];

### Variable $CONF_IND1 stores an index of the first variant in the pair
    my $CONF_IND1;

### Variable $CONF_IND2 stores an index of the second variant in the pair
    my $CONF_IND2;

    if ($CONF_IND_PAIR =~ /^(\d+)\.(\d+)$/)
    {
	$CONF_IND1 = $1;
	$CONF_IND2 = $2;
    }
    else
    {
	die "\nUnexpected format of the field supposed to contain an index of a pair of conflicting variants: $CONF_IND_PAIR\nat line $TOTAL_CONFLICTING_PAIRS_READ: $_. Stopped.";
    }

### Put an index of the pair of conflicting variants in the hash ($CONFLICTING_PAIRS_INDICES_HASH)
    $CONFLICTING_PAIRS_INDICES_HASH->{$CONF_IND_PAIR} = 1;

### Put an index of the first variant in the pair in the hash ($CONFLICTING_SINGLE_INDICES_HASH)
    $CONFLICTING_SINGLE_INDICES_HASH->{$CONF_IND1} = 1;

### Put an index of the second variant in the pair in the hash ($CONFLICTING_SINGLE_INDICES_HASH)
    $CONFLICTING_SINGLE_INDICES_HASH->{$CONF_IND2} = 1;

### Add index of the pair of conflicting variants to the array containing a list of such pairs for the first variant
    push @{$CONFLICTING_SINGLE_INDICES_TO_PAIRS->{$CONF_IND1}},$CONF_IND_PAIR;

### Add index of the pair of conflicting variants to the array containing a list of such pairs for the second variant
    push @{$CONFLICTING_SINGLE_INDICES_TO_PAIRS->{$CONF_IND2}},$CONF_IND_PAIR;

}

close(CONFIND);

my $TOTAL_CONFLICTING_PAIRS_INDICES_HASH = keys(%{$CONFLICTING_PAIRS_INDICES_HASH});
my $TOTAL_CONFLICTING_SINGLE_INDICES_HASH = keys(%{$CONFLICTING_SINGLE_INDICES_HASH});

print "\n";
print "Total number of unique indices of conflicting variant pairs:\n";
print "$TOTAL_CONFLICTING_PAIRS_INDICES_HASH\n";

print "\n";
print "Total number of unique indices of variants involved in any conflicting pair:\n";
print "$TOTAL_CONFLICTING_SINGLE_INDICES_HASH\n";
print "\n";
print "#######################################################################################################\n";


### Proceed to reading the input file with phased haplotype blocks generated with HapCUT2

print "Input file with phased haplotype blocks assembled with HapCUT2:\n$haplotype_file\n";
print "\n";

### Declare variables used in this section of the script

### Each time a new haplotype block is read, variable $total_blocks is incremented 
### When the script has finished reading the input haplotype file generated with HapCUT2, variable $total_blocks stores the total number of phased haplotype blocks present in the input file 
my $total_blocks = 0;

### Number of phased variants belonging to the currently processed haplotype block as declared in the block header line
my $block_phased_variants_parse = 0;

### Number of phased variants belonging to the currently processed haplotype block as determined from processing lines describing variants belonging to this block 
my $block_phased_variants_read = 0;

### Header line for the currently processed haplotype block
my $block_full_id;

my $blocks_lowqual_variants_hash;

### Variable $blocks_conflicting_variants_hash is a reference to a hash storing indices of blocks harboring variants involved in any conflicting pair
### It is autovivified to be a hash reference, when the first element is added
### Indices of blocks harboring variants involved in conflicting pairs are stored as hash keys
my $blocks_conflicting_variants_hash;

### Variable $blocks_conflicting_lowqual_variants_hash is a reference to a nested hash storing indices of blocks harboring variants involved in conflicting pairs and/or containing problematic variants
### It is autovivified, when the first element is added
### Any variant with mismatch or switch quality < $EXP_QUAL_VALUE is considered to be problematic
my $blocks_conflicting_lowqual_variants_hash;

### A threshold value of switch and mismatch quality to consider a variant phased with high confidence
### Any variant with mismatch or switch quality < $EXP_QUAL_VALUE is considered to be problematic (possible phasing error)
my $EXP_QUAL_VALUE = 100;

### Pruning status value for those variants that were not pruned
my $EXP_PRUN_STAT = 0;

### Variable $blocks_all_variants_lists is a reference to a nested hash containing lists of variant indices included in haplotype blocks 
### It is autovivified to be a hash reference, when the first element is added
my $blocks_all_variants_lists;

### Variable $blocks_lowqual_variants_lists is a reference to a nested hash containing per block information on problematic variants belonging to the given block 
### It is autovivified to be a hash reference, when the first element is added
my $blocks_lowqual_variants_lists;

### Total number of phased variants (included in any phased block)
my $total_phased_variants_read = 0;

### Total number of pruned variants - this number is usually very low
my $total_phased_variants_read_pruned = 0;

### Total number of phased variants involved in any conflicting pair
my $total_phased_variants_read_conflicts = 0;

### Variable $block_numbers_to_ids_map is a reference to a hash mapping haplotype block index numbers to haplotype block headers 
### It is autovivified to be a hash reference, when the first element is added
my $block_numbers_to_ids_map;

### Variable $block_numbers_to_split_ids_map is a reference to a hash mapping chunks of haplotype blocks resulting from the split of original blocks to the corresponding block chunk headers 
### It is autovivified to be a hash reference, when the first element is added
my $block_numbers_to_split_ids_map;

## Variable $map_variants_indices_to_genomic_coordinates is a reference to a nested hash mapping VCF indices of variants included in a particular haplotype block to their genomic positions
### It is autovivified to be a hash reference, when the first element is added
my $map_variants_indices_to_genomic_coordinates;

### Variable $blocks_phased_data_hash is a reference to a nested hash containing information on variants included in particular blocks
### It is autovivified to be a hash reference, when the first element is added
my $blocks_phased_data_hash;

### Variable $split_blocks_phased_data_hash is a reference to a nested hash which will be used to store information on variants included in the chunks of original phased blocks retained after the split of original blocks at problematic positions 
### It is autovivified to be a hash reference, when the first element is added
my $split_blocks_phased_data_hash;

#######################################################################################################
### Read the input file with haplotype blocks generated by HapCUT2
### Information on the format of files with haplotype blocks produced with HapCUT2 can be found at the following link:
### https://github.com/vibansal/HapCUT2

open (HPF,"$haplotype_file") or die $!;

while (<HPF>)
{
    chomp();

### Read a block header - marks start of a new block
    if (/^BLOCK\:/)
    {
	$block_full_id = $_;

### Increment the total number of haplotype blocks read from the haplotype file
### While reading the haplotype file, $total_blocks corresponds to the index number of the current haplotype block 
	$total_blocks++;

	# print "\nStart of the block:\n";
	# print "$_\n";

### Check format of the block header line
	unless (/^BLOCK\:\s+offset\:\s+(\d+)\s+len\:\s+(\d+)\s+phased\:\s+(\d+)\s+SPAN\:\s+(\d+)\s+fragments\s+(\d+)$/)
	{
	    die "\nNew haplotype block header line does not correspond to the expected format: $_. Stopped.";
	}

### Extract total number of phased variants in the block
	$block_phased_variants_parse = $3;

### Check that the hash storing header lines for haplotype blocks does not yet have an element with a key corresponding to the index number of this haplotype block
### The corresponding hash is accessed via a hash reference ($block_numbers_to_ids_map)
### Store the header line for the new block: this is done by adding a new element to the hash  
### Assign a new element to the hash: the hash key - index number of the new haplotype block, the hash value - its header line

	unless (exists($block_numbers_to_ids_map->{$total_blocks}))
	{ $block_numbers_to_ids_map->{$total_blocks} = $block_full_id; }
	else       
	{ die "\nAn element with the key $total_blocks has already been added to the hash (hash reference \$block_numbers_to_ids_map). Stopped.";}

    }

### Process lines marking the end of a haplotype block 
    elsif (/^\*{8}\s+$/)
    {
	unless (defined($block_phased_variants_parse))
	{
	    die "\nNumber of phased variants supposed to be declared in the block header line was not determined.\nCurrent block header line: $block_full_id\nCheck the format of block header lines. Stopped.";
	}

### Check that the number of phased variants in the block determined by processing lines describing individual variants is equal to the number declared in the block header line
	unless ($block_phased_variants_parse == $block_phased_variants_read)
	{
	    die "\nNumbers of phased variants declared in the block header line and determined by processing lines describing individual variants do not match: $block_phased_variants_parse vs $block_phased_variants_read.\nCurrent block header line: $block_full_id\nStopped.";
	}

### Before processing the next haplotype block, undefine the values of variables $block_phased_variants_parse, $block_full_id
	undef $block_phased_variants_parse;
	undef $block_full_id;

### Before processing the next haplotype block, set the value of the variable $block_phased_variants_read to 0
	$block_phased_variants_read = 0;
    }

### Process lines describing individual phased variants included in the current haplotype block
    else 
    {
        unless (defined($block_full_id))
	{
	    die "\nProcessing lines describing phased variants, while the block header line is not set.\nThis should not happen. Could not determine to which block a variant is assigned.\nCheck the format of the input file with haplotype blocks. Stopped.";
	}

### Increment the number of phased variants belonging to the currently processed haplotype block
	$block_phased_variants_read++;

### Increment the overall number of phased variants (included in any phased block)
	$total_phased_variants_read++;

	my @PHASED_SPLIT_FIELDS = split/\s+/;

### Get the number of elements (fields) in the currently processed line
	my $dim_split_fields = scalar(@PHASED_SPLIT_FIELDS);

### Check that the number of elements (fields) in the currently processed line is equal to the expected number (11)
### Expected number of fields is as specified in the description of the HapCUT2 output format at: https://github.com/vibansal/HapCUT2 
	unless ($dim_split_fields == 11)
	{
	    die "\nUnexpected number of fields in the line describing phased variants: $dim_split_fields vs 11 expected.\nCurrent line: $_\nCurrent block header line: $block_full_id\nCheck the format of the input file with haplotype blocks. Stopped.";
	}

### Field 1 is index of the variant in the VCF file supplied to HapCUT2
### Index of a variant is just its index number in the VCF file 
### Index of a variant is not its chromosomal position!!!
	my $PHASED_VAR_INDEX = $PHASED_SPLIT_FIELDS[0];

### Field 4 is chromosome/contig
	my $PHASED_VAR_CTG = $PHASED_SPLIT_FIELDS[3];
### Field 5 is the variant position
	my $PHASED_VAR_POS = $PHASED_SPLIT_FIELDS[4];

### Remember that the current variant with the index $PHASED_VAR_INDEX is included in the haplotype block with the index $total_blocks
	$blocks_all_variants_lists->{$total_blocks}->{$PHASED_VAR_INDEX} = 1;

### Add the element corresponding to the currently processed phased variant to a nested hash containing information on variants included in particular blocks
### $blocks_phased_data_hash is a reference to the nested hash 
### $total_blocks is the index number of the current haplotype block
### $PHASED_VAR_INDEX is the VCF index of the currently processed phased variant
	unless (exists($blocks_phased_data_hash->{$total_blocks}->{$PHASED_VAR_INDEX}))
	{
	    $blocks_phased_data_hash->{$total_blocks}->{$PHASED_VAR_INDEX} = $_;
	}
	else
	{
	    die "\nVariant with the VCF index $PHASED_VAR_INDEX included in the block with the index $total_blocks has already been added to the hash (hash reference \$blocks_phased_data_hash)\nCurrent line: $_\nCurrent block header line: $block_full_id\nVariants with duplicate indices are not supposed to be present. Stopped.;"
	}

### Store information on the genomic position for the current phased variant
	$map_variants_indices_to_genomic_coordinates->{$total_blocks}->{$PHASED_VAR_INDEX}->{'CTG'} = $PHASED_VAR_CTG;
	$map_variants_indices_to_genomic_coordinates->{$total_blocks}->{$PHASED_VAR_INDEX}->{'POS'} = $PHASED_VAR_POS;

### Check whether this phased variant is found among those involved in conflicting variant pairs specified in the input file with indices of conflicting variants pairs  
	if (exists($CONFLICTING_SINGLE_INDICES_HASH->{$PHASED_VAR_INDEX}))
	{

### If the variant is involved in any conflicting pair, increment the counter $total_phased_variants_read_conflicts
	    $total_phased_variants_read_conflicts++;

### Remember that the current haplotype block (with index $total_blocks) has a variant included in a conflicting pair of variants
	    $blocks_conflicting_variants_hash->{$total_blocks} = 1;

	    $blocks_conflicting_lowqual_variants_hash->{'CONFLICTING'}->{$total_blocks} = 1;
	    $blocks_conflicting_lowqual_variants_hash->{'CONFLICTING_OR_LOWQUAL'}->{$total_blocks} = 1;

	}

### Extract pruning status for the variant
	my $PRUN_STAT = $PHASED_SPLIT_FIELDS[8];

### Extract switch quality
	my $SWITCH_QUAL = $PHASED_SPLIT_FIELDS[9];

### Extract mismatch quality
	my $MM_QUAL = $PHASED_SPLIT_FIELDS[10];

### If switch quality < 100, remember that the current haplotype block harbors a variant with a 'problematic' switch quality value
	if ($SWITCH_QUAL < $EXP_QUAL_VALUE)
	{
	    $blocks_conflicting_lowqual_variants_hash->{'SWITCH'}->{$total_blocks} = 1;
	}

### If mismatch quality < $EXP_QUAL_VALUE, remember that the current haplotype block harbors a variant with a 'problematic' mismatch quality value
	if ($MM_QUAL < $EXP_QUAL_VALUE)
	{
	    $blocks_conflicting_lowqual_variants_hash->{'MISMATCH'}->{$total_blocks} = 1;
	}

### Remember that the current haplotype block harbors a 'problematic' variant, if switch or mismatch quality < 100 or if the variant was pruned 
	if (($SWITCH_QUAL < $EXP_QUAL_VALUE) || ($MM_QUAL < $EXP_QUAL_VALUE) || ($PRUN_STAT != $EXP_PRUN_STAT))
	{

	    $blocks_lowqual_variants_hash->{$total_blocks} = 1;

	    $blocks_conflicting_lowqual_variants_hash->{'LOWQUAL'}->{$total_blocks} = 1;
	    $blocks_conflicting_lowqual_variants_hash->{'CONFLICTING_OR_LOWQUAL'}->{$total_blocks} = 1;

	    $blocks_lowqual_variants_lists->{'LOWQUAL'}->{$total_blocks}->{$PHASED_VAR_INDEX} = 1;

	    # print "Block index number $total_blocks: $block_full_id\n";
	    # print "$_: PRUN_STAT -> $PRUN_STAT; SWITCH_QUAL -> $SWITCH_QUAL; MM_QUAL -> $MM_QUAL\n";
	    # print "\n";

	}
	unless ($PRUN_STAT == $EXP_PRUN_STAT)
	{

### If the variant was pruned, increment the counter $total_phased_variants_read_pruned
	    $total_phased_variants_read_pruned++;
	}

    }

}

close(HPF);

### Total number of blocks containing problematic variants
my $total_blocks_lowqual_variants_hash = keys(%{$blocks_lowqual_variants_hash});

### Total number of blocks containing variants involved in any conflicting pair
my $total_blocks_conflicting_variants_hash = keys(%{$blocks_conflicting_variants_hash});

### Identify blocks simultaneously containing variants involved in conflicting pairs and problematic variants
### This is used only to gather statistics on the number of such blocks
my @blocks_conflicting_and_lowqual = grep { my $CONF_BLOCK = $_; 
		   if ( grep { $CONF_BLOCK == $_ } keys(%{$blocks_conflicting_lowqual_variants_hash->{'LOWQUAL'}}))
		   {
		       1;
		   }
		   else
		   {
		       0;
		   }
} keys(%{$blocks_conflicting_lowqual_variants_hash->{'CONFLICTING'}});

my $total_blocks_conflicting_and_lowqual_variants = scalar(@blocks_conflicting_and_lowqual);

%{$blocks_conflicting_lowqual_variants_hash->{'CONFLICTING_AND_LOWQUAL'}} = map { $_ => 1 } @blocks_conflicting_and_lowqual;


print "##################################################\n";
print "Statistics on the numbers of phased haplotype blocks:\n";
print "\n";
print "Total number of phased haplotype blocks present in the input haplotype file (produced with HapCUT2):\n";
print "$total_blocks\n";

print "\n";
print "Total number of haplotype blocks containing variants involved in any conflicting pair:\n";
print "$total_blocks_conflicting_variants_hash\n";

print "\n";
print "Total number of haplotype blocks containing problematic variants:\n";
print "$total_blocks_lowqual_variants_hash\n";

print "\n";
print "Total number of haplotype blocks simultaneously containing problematic variants and variants involved in any conflicting pair:\n";
print "$total_blocks_conflicting_and_lowqual_variants\n";
print "\n";

print "\n";
print "##################################################\n";
print "Statistics on the numbers of phased variants included in haplotype blocks:\n";
print "\n";
print "Total number of phased variants (included in any phased block):\n";
print "$total_phased_variants_read\n";
print "\n";
print "Total number of phased variants involved in any conflicting pair:\n";
print "$total_phased_variants_read_conflicts\n";
print "\n";
print "Total number of pruned variants:\n";
print "$total_phased_variants_read_pruned\n";
print "\n";

print "\n";
print "#######################################################################################################\n";
print "#######################################################################################################\n";
print "Further process blocks comprising problematic variants (those with mismatch or switch qualities < $EXP_QUAL_VALUE):\n";
print "\n";

#######################################################################################################
### The next block of code does the following:
### Iterates over blocks comprising problematic variants
### Skips blocks carrying variants involved in conflicts and blocks harboring more than one problematic variant
### Tries to split remaining blocks with a single problematic variant at the position of the problematic variant
### If any chunk of the block resulting from the split has at least two phased variants retained, defines a new phased block corresponding to this chunk 
#######################################################################################################

### Variable $NBLOCK_WEIRD_VARS_DISTR is a reference to a hash storing distribution of per block numbers of problematic variants 
my $NBLOCK_WEIRD_VARS_DISTR;

### Total number of looked up blocks harboring problematic variants
my $NBLOCKS_LOWQUAL_LOOKED_UP = 0; 

### Total number of looked up blocks harboring both problematic variants and variants involved in conflicts
### Such blocks as well as blocks harboring only variants involved in conflicts will be skipped entirely
my $NBLOCKS_LOWQUAL_LOOKED_UP_CONFLICTS = 0; 

### Total number of looked up blocks harboring problematic variants but free of variants involved in conflicts
### Some of these blocks will be split at a problematic variant and the resulting flanks retained for the subsequent analysis 
my $NBLOCKS_LOWQUAL_LOOKED_UP_NO_CONFLICTS = 0; 

### Among those blocks free of variants involved in conflicts, number of those carrying a single problematic variant
my $NBLOCKS_LOWQUAL_SINGLE_WEIRD = 0;

### A reference to a hash storing numbers of split blocks for which both chunks were retained, only left chunk was retained, etc.
my $NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS;

### Set the minimum number of phased variants in a block chunk for the chunk to be retained
my $MIN_VAR_KEEP_FLANK = 2;

### Iterate over blocks comprising problematic variants
WBLOCKS: foreach my $NEW_BLOCK_TO_LOOK ( sort {$a<=>$b} keys(%{$blocks_lowqual_variants_lists->{'LOWQUAL'}}) )
{

    $NBLOCKS_LOWQUAL_LOOKED_UP++;

### Skip blocks comprising variants involved in conflicts
    if (exists($blocks_conflicting_lowqual_variants_hash->{'CONFLICTING'}->{$NEW_BLOCK_TO_LOOK}))
    {
	$NBLOCKS_LOWQUAL_LOOKED_UP_CONFLICTS++;
	next WBLOCKS; 
    }
    else
    {
	$NBLOCKS_LOWQUAL_LOOKED_UP_NO_CONFLICTS++;
    }

### Retrieve a list of indices of problematic variants found in this block
    my @BLOCK_WEIRD_VARS_INDICES = (sort {$a<=>$b} keys(%{$blocks_lowqual_variants_lists->{'LOWQUAL'}->{$NEW_BLOCK_TO_LOOK}}));

### Get the total number of problematic variants found in this block 
    my $NBLOCK_WEIRD_VARS = scalar(@BLOCK_WEIRD_VARS_INDICES);

    $NBLOCK_WEIRD_VARS_DISTR->{$NBLOCK_WEIRD_VARS}++;

### Retrieve the header line for this block
    my $NEW_BLOCK_TO_LOOK_FULLID = $block_numbers_to_ids_map->{$NEW_BLOCK_TO_LOOK};

    unless ($NEW_BLOCK_TO_LOOK_FULLID =~ /^BLOCK\:\s+offset\:\s+(\d+)\s+len\:\s+(\d+)\s+phased\:\s+(\d+)\s+SPAN\:\s+(\d+)\s+fragments\s+(\d+)$/)
    {
	    die "\nHaplotype block header for the block $NEW_BLOCK_TO_LOOK does not correspond to the expected format: $NEW_BLOCK_TO_LOOK_FULLID. Stopped.";
    }

    my $NEW_BLOCK_OFFSET = $1;
    my $NEW_BLOCK_LENGTH = $2;
    my $NEW_BLOCK_PHASED = $3;
    my $NEW_BLOCK_SPAN = $4;

### Process the block with problematic variants further only if it carries a single problematic variant
    if ($NBLOCK_WEIRD_VARS == 1)
    {
	$NBLOCKS_LOWQUAL_SINGLE_WEIRD++;

### Retrieve the list of indices for all variants belonging to this block
	my @BLOCK_ALL_VARS_INDICES = (sort {$a<=>$b} keys(%{$blocks_all_variants_lists->{$NEW_BLOCK_TO_LOOK}}));

	my $MIN_PHASED_INDEX = min(@BLOCK_ALL_VARS_INDICES);
	my $MAX_PHASED_INDEX = max(@BLOCK_ALL_VARS_INDICES);

### Compute length of the block expressed as the number of encompassed variants (both phased and unphased)
	my $NEW_BLOCK_LENGTH_COMP = $MAX_PHASED_INDEX-$MIN_PHASED_INDEX+1; 

### Get the total number of phased variants belonging to this block
	my $NBLOCK_ALL_VARS = scalar(@BLOCK_ALL_VARS_INDICES);

	my $MIN_PHASED_INDEX_POS = $map_variants_indices_to_genomic_coordinates->{$NEW_BLOCK_TO_LOOK}->{$MIN_PHASED_INDEX}->{'POS'};
	my $MAX_PHASED_INDEX_POS = $map_variants_indices_to_genomic_coordinates->{$NEW_BLOCK_TO_LOOK}->{$MAX_PHASED_INDEX}->{'POS'};

### Compute length of the block in the usual way, as the number of spanned genomic positions
	my $NEW_BLOCK_SPAN_COMP = $MAX_PHASED_INDEX_POS-$MIN_PHASED_INDEX_POS; 

	# print "$NEW_BLOCK_TO_LOOK_FULLID: @BLOCK_ALL_VARS_INDICES\n";
	# print "$NEW_BLOCK_TO_LOOK_FULLID: $MIN_PHASED_INDEX\n";
	# print "$NEW_BLOCK_TO_LOOK_FULLID: $NEW_BLOCK_LENGTH_COMP vs $NEW_BLOCK_LENGTH\n";
	# print "$NEW_BLOCK_TO_LOOK_FULLID: $NBLOCK_ALL_VARS vs $NEW_BLOCK_PHASED\n";
	# print "NEW_BLOCK_TO_LOOK: $NEW_BLOCK_TO_LOOK; $NEW_BLOCK_TO_LOOK_FULLID; NBLOCK_ALL_VARS $NBLOCK_ALL_VARS\n";

	# print "$NEW_BLOCK_TO_LOOK_FULLID: MIN_PHASED_INDEX $MIN_PHASED_INDEX -> $MIN_PHASED_INDEX_POS\n";
	# print "$NEW_BLOCK_TO_LOOK_FULLID: MAX_PHASED_INDEX $MAX_PHASED_INDEX -> $MAX_PHASED_INDEX_POS\n";
	# print "$NEW_BLOCK_TO_LOOK_FULLID: NEW_BLOCK_SPAN_COMP $NEW_BLOCK_SPAN_COMP\n";
	# print "\n";

	unless ($MIN_PHASED_INDEX == $NEW_BLOCK_OFFSET)
	{
	    die "\nIndex of the first phased variant in the block $NEW_BLOCK_TO_LOOK does not match the declared offset: $MIN_PHASED_INDEX vs $NEW_BLOCK_OFFSET.\nBlock header line: $NEW_BLOCK_TO_LOOK_FULLID\nStopped.";
	}

	unless ($NEW_BLOCK_LENGTH_COMP == $NEW_BLOCK_LENGTH)
	{
	    die "\nNumber of variants (phased and unphased) encompassed by the block $NEW_BLOCK_TO_LOOK does not match the declared number of encompassed variants: $NEW_BLOCK_LENGTH_COMP vs $NEW_BLOCK_LENGTH.\nBlock header line: $NEW_BLOCK_TO_LOOK_FULLID\nStopped.";
	}

	unless ($NBLOCK_ALL_VARS == $NEW_BLOCK_PHASED)
	{
	    die "\nNumber of phased variants listed for the block $NEW_BLOCK_TO_LOOK does not match the declared number of phased variants: $NBLOCK_ALL_VARS vs $NEW_BLOCK_PHASED.\nBlock header line: $NEW_BLOCK_TO_LOOK_FULLID\nStopped.";

	}
	unless ($NEW_BLOCK_SPAN_COMP == $NEW_BLOCK_SPAN)
	{
	    die "\nNumber of positions encompassed by the block $NEW_BLOCK_TO_LOOK does not match the declared number of the encompassed positions: $NEW_BLOCK_SPAN_COMP vs $NEW_BLOCK_SPAN\nBlock header line: $NEW_BLOCK_TO_LOOK_FULLID\nStopped.";
	}

### Retrieve the index of the single problematic variant found in the block
	my $BLOCK_WEIRD_VARIANT_INDEX = $BLOCK_WEIRD_VARS_INDICES[0];

	unless (grep { $_ == $BLOCK_WEIRD_VARIANT_INDEX } @BLOCK_ALL_VARS_INDICES)
	{
	    die "\nPutative problematic variant (VCF index $BLOCK_WEIRD_VARIANT_INDEX) found in the block $NEW_BLOCK_TO_LOOK is missing from the list of variants belonging to this block.\nBlock header line: $NEW_BLOCK_TO_LOOK_FULLID\nStopped.";	   
	} 

#######################################################################################################
### Split the block at the problematic variant

### Retrieve variants included in the block and located upstream of the problematic variant - these variants constitute the left flank
	my @VARIANT_INDICES_LEFT_FLANK = grep { $_ < $BLOCK_WEIRD_VARIANT_INDEX } @BLOCK_ALL_VARS_INDICES;

### Retrieve variants included in the block and located downstream of the problematic variant - these variants constitute the right flank
	my @VARIANT_INDICES_RIGHT_FLANK = grep { $_ > $BLOCK_WEIRD_VARIANT_INDEX } @BLOCK_ALL_VARS_INDICES;

	my $DIM_LEFT_FLANK = scalar(@VARIANT_INDICES_LEFT_FLANK);
	my $DIM_RIGHT_FLANK = scalar(@VARIANT_INDICES_RIGHT_FLANK);


### Check whether the flanks resulting from the split of this phased block carry sufficient numbers (>=2) of phased variants to be considered as separate blocks
### If any flank has a sufficient number of phased variants, define a new block corresponding to this flank

	my $SPLIT_INDICES_REF_HASH;

	if (($DIM_LEFT_FLANK >= $MIN_VAR_KEEP_FLANK) && ($DIM_RIGHT_FLANK >= $MIN_VAR_KEEP_FLANK))
	{
	    $NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS->{'BOTH_FLANKS'}++;

	    $SPLIT_INDICES_REF_HASH->{'LEFT_FLANK'} = \@VARIANT_INDICES_LEFT_FLANK;
	    $SPLIT_INDICES_REF_HASH->{'RIGHT_FLANK'} = \@VARIANT_INDICES_RIGHT_FLANK;

	}
	elsif ($DIM_LEFT_FLANK >= $MIN_VAR_KEEP_FLANK)
	{
	    $NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS->{'LEFT_FLANK'}++;

	    $SPLIT_INDICES_REF_HASH->{'LEFT_FLANK'} = \@VARIANT_INDICES_LEFT_FLANK;

	}
	elsif ($DIM_RIGHT_FLANK >= $MIN_VAR_KEEP_FLANK)
	{
	    $NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS->{'RIGHT_FLANK'}++;

	    $SPLIT_INDICES_REF_HASH->{'RIGHT_FLANK'} = \@VARIANT_INDICES_RIGHT_FLANK;
	}
	else
	{
	    $NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS->{'SHORT_FLANKS'}++;
	    next WBLOCKS;
	}

	$NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS->{'ANY_FLANK'}++;

### If any of the flanks is to be considered as a separate block, construct the corresponding header line and put information on this block to the hash
	foreach my $SPLIT_PART1 ( sort {$a cmp $b} keys(%{$SPLIT_INDICES_REF_HASH}) )
	{

### Header line for the new block resulting from the split is constructed using subroutine 'define_split_blocks'
	    my $NEW_SPLIT_PART_ID = define_split_blocks($SPLIT_INDICES_REF_HASH->{$SPLIT_PART1}, $map_variants_indices_to_genomic_coordinates, $NEW_BLOCK_TO_LOOK);

### Put the header line constructed for the block chunk to the hash mapping it to the original block 
	    $block_numbers_to_split_ids_map->{$NEW_BLOCK_TO_LOOK}->{$SPLIT_PART1} = $NEW_SPLIT_PART_ID;

	    my @NEW_SPLIT_PART_INDICES = @{$SPLIT_INDICES_REF_HASH->{$SPLIT_PART1}};

### Iterate over all variants included in this chunk of the original block
	    foreach my $NEW_SPLIT_PART_INDEX (@NEW_SPLIT_PART_INDICES)
	    {
		unless (exists($blocks_phased_data_hash->{$NEW_BLOCK_TO_LOOK}->{$NEW_SPLIT_PART_INDEX}))
		{
		    die "\nVariant with the index $NEW_SPLIT_PART_INDEX belonging to the newly defined block chunk is missing from the list of variants included in the original block (block index $NEW_BLOCK_TO_LOOK).\nBlock chunk header line: $NEW_SPLIT_PART_ID\nOriginal block header line: $NEW_BLOCK_TO_LOOK_FULLID. Stopped.";
		}

### Retrieve the string with phasing information for this variant as was specified in the original block
		my $NEW_SPLIT_PART_INDEX_STRING = $blocks_phased_data_hash->{$NEW_BLOCK_TO_LOOK}->{$NEW_SPLIT_PART_INDEX};

### Put this string into to the hash containing phasing information for variants included in individual block chunks retained after the split of original blocks 
		$split_blocks_phased_data_hash->{$NEW_BLOCK_TO_LOOK}->{$SPLIT_PART1}->{$NEW_SPLIT_PART_INDEX} = $NEW_SPLIT_PART_INDEX_STRING;
	    }

	}

    }

### Skip a block if it harbors more than one problematic variant
    else
    {
	next WBLOCKS;
    }

}

print "Among $NBLOCKS_LOWQUAL_LOOKED_UP phased blocks harboring one or more problematic variants:\n";
print "$NBLOCKS_LOWQUAL_LOOKED_UP_CONFLICTS harbor variants involved in any conflicting pair (these were skipped)\n";
print "$NBLOCKS_LOWQUAL_LOOKED_UP_NO_CONFLICTS are free of variants involved in conflicting pairs (these were processed further)\n";

print "\n";
print "Per block distribution of numbers of problematic variants among blocks harboring one or more problematic variants\n(only $NBLOCKS_LOWQUAL_LOOKED_UP_NO_CONFLICTS blocks free of conflicting variants are considered):\n";
print "Number of problematic variants - number of blocks\n", Dump $NBLOCK_WEIRD_VARS_DISTR;

print "\n";
print "##################################################\n";
print "$NBLOCKS_LOWQUAL_SINGLE_WEIRD blocks harboring exactly one problematic variant were subjected to further analysis.\n";
print "These blocks were split at the position of the problematic variant.\n";
print "If a chunk resulting from the split had at least two phased variants, it was retained and treated as a separated block.\n";

print "\n";
print "For $NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS->{'ANY_FLANK'} blocks, at least one of the two chunks was retained.\n";
print "For $NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS->{'SHORT_FLANKS'} blocks, both chunks were discarded as neither of the two chunks had a sufficient number of variants retained.\n";

### The next commented line allows to print out distribution of blocks (both flanks retained, only left flank retained, etc.)
# print "\nWhich flank retained: number of split blocks\n", Dump $NBLOCKS_LOWQUAL_KEEP_FLANKS_COUNTS;


#######################################################################################################
### Output filtered haplotype blocks

print "\n";
print "#######################################################################################################\n";
print "#######################################################################################################\n";
print "\n";
print "Report filtered haplotype blocks:\n";

### String to use as a haplotype block closing line according to the HapCUT2 output format
my $BLOCK_CLOSING_LINE = '********';

#######################################################################################################
### 1) Report an updated list of phased blocks subjected to filtering based only on the presence of variants involved in conflicts (no filtering for 'problematic' variants applied)
### The resulting file contains those original phased blocks that are free of variants involved in conflicts

my $NBLOCKS_ANALYZED_TOTAL1 = 0;
my $NBLOCKS_ANALYZED_WO_CONFLICTS = 0;
my $NBLOCKS_ANALYZED_W_CONFLICTS = 0;

open (HPCFILT,">$haplotype_file_conflicts_filtered") or die $!;

PBLOCKS1: foreach my $ORIGINAL_BLOCK1 (sort {$a<=>$b} keys(%{$blocks_phased_data_hash}))
{
    $NBLOCKS_ANALYZED_TOTAL1++;

    my $ORIGINAL_BLOCK_FULLID1 = $block_numbers_to_ids_map->{$ORIGINAL_BLOCK1};

    unless (exists($blocks_conflicting_lowqual_variants_hash->{'CONFLICTING'}->{$ORIGINAL_BLOCK1}))
    {
	$NBLOCKS_ANALYZED_WO_CONFLICTS++;

	print HPCFILT "$ORIGINAL_BLOCK_FULLID1\n";

	foreach my $FINE_BLOCK_INDEX1 (sort {$a<=>$b} keys(%{$blocks_phased_data_hash->{$ORIGINAL_BLOCK1}}))
	{
	    my $FINE_BLOCK_INDEX_LINE1 = $blocks_phased_data_hash->{$ORIGINAL_BLOCK1}->{$FINE_BLOCK_INDEX1};
	    print HPCFILT "$FINE_BLOCK_INDEX_LINE1\n";
	}

	print HPCFILT "$BLOCK_CLOSING_LINE\n";

    }
    else
    {
	$NBLOCKS_ANALYZED_W_CONFLICTS++;
    }

}

close(HPCFILT);

print "\n";
print "#######################################################################################################\n";
print "1. Haplotype blocks filtered based on the list of variants involved in conflicting variant pairs were written to the file:\n";
print "$haplotype_file_conflicts_filtered\n";

print "\n";
print "This output file contains those original haplotype blocks that were free from variants involved in conflicting variant pairs (no further filtering based on HapCUT2 phasing qualities was applied).\n";

print "\n";
print "Total number of original phased haplotype blocks present in the input haplotype file:\n";
print "$total_blocks\n";

print "\n";
print "Total number of filtered out haplotype blocks (containing variants involved in any conflicting pair):\n";
print "$total_blocks_conflicting_variants_hash\n";

print "\n";
print "Total number of haplotype blocks in the resulting output file:\n";
system("grep -c '^BLOCK:' $haplotype_file_conflicts_filtered");

print "#######################################################################################################\n";

#######################################################################################################
### 2) Report an updated list of phased blocks subjected to filtering based both based on the supplied list of conflicting variant indices and based on probabilities of phasing errors computed with HapCUT2
### The resulting file contains:
### a) original intact phased blocks that are free both of variants involved in conflicts and problematic variants
### b) chunks of original blocks resulting from the split of an original block at the position of a problematic variant
### only those blocks that were devoid of variants involved in conflicts and carried a single problematic variant were split

### Total number of processed phased haplotype blocks (those present in the input file) 
my $NBLOCKS_ANALYZED_TOTAL = 0;

### Total number of processed blocks free of variants involved in conflicts and problematic variants - they are printed intact to the output file
my $NBLOCKS_ANALYZED_FINE = 0;

### Total number of processed blocks carrying variants involved in conflicts and/or problematic variants
my $NBLOCKS_ANALYZED_FAULTY = 0;

### Total number of processed blocks carrying variants involved in conflicts
my $NBLOCKS_ANALYZED_CONFLICTS = 0;

## Total number of processed blocks devoid of variants involved in conflicts, but possessing problematic variants
my $NBLOCKS_ANALYZED_LOWQUAL = 0;

## Total number of processed blocks with a problematic variant, for which at least one chunk was retained after the split
my $NBLOCKS_ANALYZED_LOWQUAL_SPLIT = 0;

### Total number of chunks of original phased blocks to output (resulting from the split of original blocks at a problematic position)
my $NSPLIT_PARTS_TOTAL = 0;

### Total number of haplotype blocks to output - includes both intact original phased blocks and individual chunks of original blocks resulting from the split of original blocks at a problematic position
my $NVIRTUAL_BLOCKS_TO_REPORT = 0; 

open (HPFILT,">$haplotype_file_filtered") or die $!;

PBLOCKS: foreach my $ORIGINAL_BLOCK (sort {$a<=>$b} keys(%{$blocks_phased_data_hash}))
{
    $NBLOCKS_ANALYZED_TOTAL++;

    my $ORIGINAL_BLOCK_FULLID = $block_numbers_to_ids_map->{$ORIGINAL_BLOCK};

### Process blocks with variants involved in conflicts and/or problematic variants
    if (exists($blocks_conflicting_lowqual_variants_hash->{'CONFLICTING_OR_LOWQUAL'}->{$ORIGINAL_BLOCK}))
    {
	$NBLOCKS_ANALYZED_FAULTY++;

### Skip blocks comprising variants involved in conflicts - these are not included in the output file
	if (exists($blocks_conflicting_lowqual_variants_hash->{'CONFLICTING'}->{$ORIGINAL_BLOCK}))
	{
	    $NBLOCKS_ANALYZED_CONFLICTS++;
	    next PBLOCKS;
	}
### Process blocks devoid of variants involved in conflicts, but possessing problematic variants
	elsif (exists($blocks_conflicting_lowqual_variants_hash->{'LOWQUAL'}->{$ORIGINAL_BLOCK}))
	{
	    $NBLOCKS_ANALYZED_LOWQUAL++;

### Check whether this block was subjected to splitting and at least one resulting chunk was retained 
	    if (exists($split_blocks_phased_data_hash->{$ORIGINAL_BLOCK}))
	    {
    		$NBLOCKS_ANALYZED_LOWQUAL_SPLIT++;

		my @SPLIT_PARTS = sort { $a cmp $b } keys(%{$split_blocks_phased_data_hash->{$ORIGINAL_BLOCK}});

### Number of chunks retained for this block after split - can be 1 (one of the two flanks retained) or 2 (both flanks retained)
		my $NSPLIT_PARTS = scalar(@SPLIT_PARTS);
		$NSPLIT_PARTS_TOTAL += $NSPLIT_PARTS;

		$NVIRTUAL_BLOCKS_TO_REPORT += $NSPLIT_PARTS;

### Print each retained chunk of this split block to the output file 
		foreach my $SPLIT_PART (@SPLIT_PARTS)
		{
		    unless (exists($block_numbers_to_split_ids_map->{$ORIGINAL_BLOCK}->{$SPLIT_PART}))
		    {
			die "\nNo header line was constructed and specified for chunk $SPLIT_PART of the original block $ORIGINAL_BLOCK.\nOriginal block header line: $ORIGINAL_BLOCK_FULLID\nStopped.";
		    }
	
### Retrieve a header line retained for this chunk 
		    my $SPLIT_PART_FULLID = $block_numbers_to_split_ids_map->{$ORIGINAL_BLOCK}->{$SPLIT_PART};

		    print HPFILT "$SPLIT_PART_FULLID\n";

		    foreach my $SPLIT_BLOCK_PART_INDEX (sort {$a<=>$b} keys(%{$split_blocks_phased_data_hash->{$ORIGINAL_BLOCK}->{$SPLIT_PART}}))
		    {
		    	my $SPLIT_BLOCK_PART_INDEX_LINE = $split_blocks_phased_data_hash->{$ORIGINAL_BLOCK}->{$SPLIT_PART}->{$SPLIT_BLOCK_PART_INDEX};
		    	print HPFILT "$SPLIT_BLOCK_PART_INDEX_LINE\n";

		    }

		    print HPFILT "$BLOCK_CLOSING_LINE\n";

		}

	    }

	}

    }
### Print blocks devoid both of variants involved in conflicting pairs and problematic variants
    else
    {
	$NBLOCKS_ANALYZED_FINE++;
	$NVIRTUAL_BLOCKS_TO_REPORT++; 

	print HPFILT "$ORIGINAL_BLOCK_FULLID\n";

	foreach my $FINE_BLOCK_INDEX (sort {$a<=>$b} keys(%{$blocks_phased_data_hash->{$ORIGINAL_BLOCK}}))
	{
	    my $FINE_BLOCK_INDEX_LINE = $blocks_phased_data_hash->{$ORIGINAL_BLOCK}->{$FINE_BLOCK_INDEX};
	    print HPFILT "$FINE_BLOCK_INDEX_LINE\n";
	}

	print HPFILT "$BLOCK_CLOSING_LINE\n";
   }

}

close(HPFILT);

print "\n";
print "#######################################################################################################\n";
print "2. Haplotype blocks subjected to filtering based both on the supplied list of conflicting variant indices and probabilities of phasing errors computed by HapCUT2 were written to the file:\n";
print "$haplotype_file_filtered\n";

print "\n";
print "This output file contains:\n"; 
print "a) original intact haplotype blocks that were free both from variants involved in conflicts and problematic variants\n";
print "b) chunks of original haplotype blocks resulting from the split of an original block at the position of a problematic variant*\n";

print "\n*Only those blocks that were devoid of variants involved in conflicts and carried a single problematic variant were subjected to splitting\n";

print "\n";
print "Total number of original phased haplotype blocks present in the input haplotype file:\n";
print "$total_blocks\n";

print "\n";
print "Total number of haplotype blocks in the resulting output file:\n";
system("grep -c '^BLOCK:' $haplotype_file_filtered");
print "These include:\n";
print "$NBLOCKS_ANALYZED_FINE original intact haplotype blocks\n";
print "$NSPLIT_PARTS_TOTAL haplotype blocks corresponding to chunks of original blocks**\n";

print "\n**These chunks result from splitting of $NBLOCKS_ANALYZED_LOWQUAL_SPLIT original phased haplotype blocks\n";


print "#######################################################################################################\n";
print "#######################################################################################################\n";
print "\nFinished\n";
system("date");


### Subroutines

sub define_split_blocks
{
    my ($SPLIT_INDICES, $MAP_INDICES_TO_GENOMIC_COORDINATES, $PARENT_BLOCK_ID) = @_;

### Dereference a reference to an array containing indices of variants belonging to this chunk of the split block
    my @SPLIT_VARS_INDICES = @{$SPLIT_INDICES};

    my $MIN_SPLIT_PHASED_INDEX = min(@SPLIT_VARS_INDICES);
    my $MAX_SPLIT_PHASED_INDEX = max(@SPLIT_VARS_INDICES);

    my $SPLIT_BLOCK_OFFSET = $MIN_SPLIT_PHASED_INDEX; ###INDEX OF THE BLOCK START

### Compute length of this block chunk expressed as the number of encompassed variants (both phased and unphased)
    my $SPLIT_BLOCK_LENGTH = $MAX_SPLIT_PHASED_INDEX-$MIN_SPLIT_PHASED_INDEX+1; 

### Get the total number of phased variants belonging to this block chunk
    my $SPLIT_BLOCK_PHASED = scalar(@SPLIT_VARS_INDICES);

### Retrieve genomic positions for the first and last phased variants harbored by this block chunk
    my $MIN_SPLIT_PHASED_INDEX_POS = $MAP_INDICES_TO_GENOMIC_COORDINATES->{$PARENT_BLOCK_ID}->{$MIN_SPLIT_PHASED_INDEX}->{'POS'};
    my $MAX_SPLIT_PHASED_INDEX_POS = $MAP_INDICES_TO_GENOMIC_COORDINATES->{$PARENT_BLOCK_ID}->{$MAX_SPLIT_PHASED_INDEX}->{'POS'};

### Compute length of the block chunk in the usual way, as the number of spanned genomic positions
    my $SPLIT_BLOCK_SPAN = $MAX_SPLIT_PHASED_INDEX_POS-$MIN_SPLIT_PHASED_INDEX_POS; 

### Set number of fragments in the block chunk to 0, as information on the number of fragments cannot be determined from a haplotype file alone
### For compatibiliy with HapCutToVcf utility, this field cannot be set to 'NA'
    my $SPLIT_BLOCK_FRAGMENTS = '0';

### Construct a string corresponding to the header line for the newly defined block
    my $SPLIT_BLOCK_ID_TO_STRING = join(' ', 'BLOCK: offset:', $SPLIT_BLOCK_OFFSET, 'len:', $SPLIT_BLOCK_LENGTH, 'phased:', $SPLIT_BLOCK_PHASED, 'SPAN:', $SPLIT_BLOCK_SPAN, 'fragments', $SPLIT_BLOCK_FRAGMENTS);

    return($SPLIT_BLOCK_ID_TO_STRING);

}

