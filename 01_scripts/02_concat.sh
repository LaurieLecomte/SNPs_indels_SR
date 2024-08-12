#!/bin/bash

# Concatenate VCFs for all chromosomes together
# Run this script from SNPs_indels_SR directory using this command :

# srun -c 2 --mem=10G -p ibis_medium --time=7-00:00:00 -J 02_concat -o log/02_merge_%j.log 01_scripts/02_concat.sh &

# VARIABLES
CHR_LIST="02_infos/chr.txt"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
MERGED_VCF="$MERGED_DIR/SNPs_indels.vcf.gz"
VCF_LIST="02_infos/VCF_list.txt"

CPU=2

# LOAD REQUIRED MODULES
module load bcftools/1.16


# 1. Make a list of VCFs to concatenate
ls -1 $CALLS_DIR/*.vcf.gz > $VCF_LIST


# 2. Concatenate files 
bcftools concat -f $VCF_LIST -Oz --threads $CPU > $MERGED_VCF