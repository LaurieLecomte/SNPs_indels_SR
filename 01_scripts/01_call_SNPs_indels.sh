#!/bin/bash

# Call SNPs and indels in all samples, parallelized by chromosome
# Run this script from SNPs_indels_SR directory using this command :
# parallel -a 02_infos/chr.txt -k -j 10 srun -c 2 --mem=10G -p ibis_medium --time=7-00:00:00 -J 01_call_SNPs_indels_{} -o log/01_call_SNPs_indels_{}_%j.log 01_scripts/01_call_SNPs_indels.sh {} &


# VARIABLES
REF="03_genome/genome.fasta"     # reference genome
BAM_LIST="02_infos/bam_list.txt" # list of bam files
CHR=$1                           # region or chromosome for this job 
CALLS_DIR="05_calls"              # results directory
CPU=2                            # number of threads to use per job

# LOAD REQUIRED MODULES
module load bcftools

# 1. Call SNPs only : required -I option
bcftools mpileup -Ou -f $REF -r $CHR -b $BAM_LIST -a AD,DP,SP,ADF,ADR -q 5 --threads $CPU | bcftools call -a GP,GQ -mv -Oz --threads $CPU > $OUT_DIR/"$CHR".vcf.gz 
