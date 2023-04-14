#!/bin/bash

# Filter SNPs 
# VERY IMPORTANT : deactivate conda env first because setGT will not work if older bcftools version is loaded by conda by default !

# Run this script from SNPs_indels_SR directory using this command :
# srun -c 4 --mem=10G -p medium --time=7-00:00:00 -J 03_filter_samples -o log/03_filter_samples_%j.log 01_scripts/03_filter_samples.sh &

# VARIABLES
CHR_LIST="02_infos/chr.txt"
BAM_LIST="02_infos/bam_list.txt"

CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
MERGED_VCF="$MERGED_DIR/SNPs_indels.vcf.gz"
VCF_LIST="02_infos/VCF_list.txt"
FILT_DIR="07_filtered"

CPU=4

REGIONS_EX="02_infos/excl_chrs.txt"

## Filtering criteria : 
MIN_GQ=5
MAX_GQ=127 # deduct max allowed GQ from VCF or plot, e.g. bcftools view -s '13066A' | bcftools query -f "[%GQ\n]" | sort -n | uniq | tail -n1

MIN_DP=4
MAX_DP=80 # arbitrary threshold of 5 * anticipated SR sequencing coverage, 5 * 16 


# LOAD REQUIRED MODULES
module load bcftools/1.15


# 1. Remove non biallelic sites, assign missing genotype where DP and GQ are too low or too extreme and remove unwanted contigs from header 
bcftools view -m2 -M2 --threads $CPU $MERGED_VCF | bcftools +setGT --threads $CPU -- -t q -n . -e "FORMAT/DP >= $MIN_DP & FORMAT/DP < $MAX_DP & FORMAT/GQ >= $MIN_GQ & FORMAT/GQ < $MAX_GQ" | bcftools sort | grep -vFf $REGIONS_EX > $FILT_DIR/SNPs_indels_DP"$MIN_DP"_GQ"$MIN_GQ".vcf

# 2. Compress and index
bgzip $FILT_DIR/SNPs_indels_DP"$MIN_DP"_GQ"$MIN_GQ".vcf -f --threads $CPU
tabix -p vcf $FILT_DIR/SNPs_indels_DP"$MIN_DP"_GQ"$MIN_GQ".vcf.gz -f


