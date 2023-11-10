#!/bin/bash

# manitou
# srun -c 2 -p small --time=1-00:00:00 --mem=50G -J SNPs_compare_FMISS_filters -o log/SNPs_compare_FMISS_filters_%j.log /bin/sh ./01_scripts/utils/compare_FMISS_filters_SNPs.sh &

# VARIABLES

# First 2 files were produced by the "utils/extract_unfiltered.sh" script
RAW_SNP="/project/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/SNPs.vcf.gz"
#RAW_INDEL="/project/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/indels.vcf.gz"

FILT_DIR="07_filtered/SNPs"

REGIONS_EX="02_infos/excl_chrs.txt"

MIN_GQ=5
MAX_GQ=256 # deduct max allowed GQ from VCF, e.g. bcftools query -f "[%GQ\n]" $CALLS_DIR/raw/"$SAMPLE".vcf | sort -n | uniq | tail -n1

MIN_DP=4
MAX_DP=80 # arbitrary threshold of 5 * anticipated SR sequencing coverage, 5 * 16 

MIN_MAF=0.05
MAX_MAF=0.95

CPU=2

#MAX_MISS=0.5

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load htslib/1.13

# 0. Add tags
bcftools +fill-tags --threads $CPU $RAW_SNP -Oz -- -t AC,AC_Hom,AC_Het,AC_Hemi,AF,AN,NS,ExcHet,HWE,MAF,F_MISSING > $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_tagged.vcf.gz

# 1. Filter on MAF and proportion of missing genotypes, with F_MISS = 0.75, 0.625, 0.125, 0.25, 0.375 
MAX_MISS=0.75
bcftools view --max-alleles 2 $RAW_SNP | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
echo "with $MAX_MISS : $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"

MAX_MISS=0.625
bcftools view --max-alleles 2 $RAW_SNP | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
echo "with $MAX_MISS : $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"


MAX_MISS=0.375
bcftools view --max-alleles 2 $RAW_SNP | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
echo "with $MAX_MISS : $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"


MAX_MISS=0.375
bcftools view --max-alleles 2 $RAW_SNP | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
echo "with $MAX_MISS : $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"


MAX_MISS=0.25
bcftools view --max-alleles 2 $RAW_SNP | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
echo "with $MAX_MISS : $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"

MAX_MISS=0.125
bcftools view --max-alleles 2 $RAW_SNP | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
echo "with $MAX_MISS : $FILT_DIR/"$(basename -s .vcf.gz $RAW_SNP)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"

