#!/usr/bin/env bash

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

INPUT_VCF_STATS=$1

if [ ! -f $INPUT_VCF_STATS ] || [ -z $INPUT_VCF_STATS ]
then 
    echo
    echo "Input file ($INPUT_VCF_STATS) does not exist or was not specified. Stopped."
    exit 1
fi

echo '#######################################################################################################'
echo 'Run compute_Fis.sh'
echo '#######################################################################################################'
echo
echo 'This script computes Fis values for individual biallelic sites based on observed and expected numbers of heterozygous genotypes.' 
echo
echo 'Fis is computed as 1-Ho/He, where Ho and He stand for the observed and expected numbers of heterozygous genotypes respectively.'
echo 'Fractional expected genotype counts are rounded to the nearest integer number.'
echo
echo 'The script takes as input a file with basic population genetics statistics produced with the populations program (http://catchenlab.life.illinois.edu/stacks/comp/populations.php) from a multi-sample VCF file.'
echo
echo "The script produces a single output file with the extension '.Fis.txt' containing computed Fis values and other related data."
echo
echo 'For further information, see README.md'
echo
echo '#######################################################################################################'
echo 'Started'
date
echo '#######################################################################################################'
echo
echo 'Input file:'
echo "$INPUT_VCF_STATS"

echo
echo 'Number of records in the input file:'
grep -cv '^#' $INPUT_VCF_STATS

#######################################################################################################
### Specify the name for the output file
declare OUTPUT_FIS_FILE

PAT=".p.sumstats.tsv$"

### Check whether the input file has the usual '.p.sumstats.tsv' extension produced by 'populations' 
### If so, create a name for the output file by replacing the '.p.sumstats.tsv' extension with the '.Fis.txt' extension
### Otherwise, construct a name for the output file by appending the '.Fis.txt' extension to the name of the input file   
if [[ $INPUT_VCF_STATS =~ $PAT ]]
then
    OUTPUT_FIS_FILE=$( echo $INPUT_VCF_STATS | sed 's/\.p\.sumstats\.tsv$/.Fis.txt/' )
else
    OUTPUT_FIS_FILE=$INPUT_VCF_STATS.Fis.txt
fi

#######################################################################################################
### Create a temporary file containing only required columns (columns 1-13) from the output of the 'populations' program

INPUT_VCF_STATS_TMP="$INPUT_VCF_STATS.tmp"

head -1 $INPUT_VCF_STATS > $INPUT_VCF_STATS_TMP
tail -n +2 $INPUT_VCF_STATS | cut -f1-13 >> $INPUT_VCF_STATS_TMP

#######################################################################################################
### Construct a header for the output file

HEADER=$(head -2 $INPUT_VCF_STATS_TMP | tail -1 | sed 's/ /\./g' | sed 's/^#\./#/')
HEADER=$(printf "$HEADER\tN.Obs.Het\tN.Exp.Het\tN.Exp.Het.rounded\tRatio_N.Obs.Het_to_N.Exp.Het.rounded\tFis")

echo "$HEADER" > $OUTPUT_FIS_FILE

#######################################################################################################
### For each variant present in the input file, compute:
### 1) Observed number of heterozygous genotypes -> column 14
### 2) Expected number of heterozygous genotypes -> column 15
### 3) Expected number of heterozygous genotypes rounded to the nearest integer number -> column 16
### 4) Ho/He, where Ho and He stand for the observed and expected numbers of heterozygous genotypes respectively. Rounded expected numbers of heterozygous genotypes are used. -> column 17
### 5) Fis computed as 1 - Ho/He. Rounded expected numbers of heterozygous genotypes are used. -> column 18
### Computed values are appended to the end of the current string.
### The updated string is printed to the output file.

grep -v '^#' $INPUT_VCF_STATS_TMP | awk 'BEGIN { OFS="\t"; } { $14=$8*$10; $15=$8*$12; $16=sprintf("%.0f",$15); $17=$14/$16; $18=(1-($14/$16)); print; }' >> $OUTPUT_FIS_FILE

echo
echo '#######################################################################################################'
echo 'Computed Fis values were written to the file:'
echo "$OUTPUT_FIS_FILE"

echo
echo 'Number of records in the output file:'
grep -cv '^#' $OUTPUT_FIS_FILE

### Delete the temporary file
rm $INPUT_VCF_STATS_TMP

echo
echo '#######################################################################################################'
echo 'Finished'
date
