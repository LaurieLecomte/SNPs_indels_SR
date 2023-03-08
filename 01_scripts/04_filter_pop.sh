#!/bin/bash

# Filter SNPs set on MAF and proportion of missing data 
# Run this script from SNPs_indels_SR directory using this command :

# manitou
# srun -c 4 --mem=10G -p medium --time=7-00:00:00 -J 04_filter_pop -o log/04_filter_pop_%j.log 01_scripts/04_filter_pop.sh &

# valeria
# srun -c 4 --mem=10G -p ibis_medium --time=7-00:00:00 -J 04_filter_pop -o log/04_filter_pop_%j.log 01_scripts/04_filter_pop.sh &

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

MIN_MAF=0.05
MAX_MISS=0.5

# LOAD REQUIRED MODULES
module load bcftools/1.13



# 1. Add tags
bcftools +fill-tags $FILT_DIR/SNPs_indels_DP"$MIN_DP"_GQ"$MIN_GQ".vcf -Oz -- -t all > $FILT_DIR/SNPs_indels_DP"$MIN_DP"_GQ"$MIN_GQ".tagged.vcf.gz

# 2. Filter on MAF and proportion of missing genotypes
bcftools view --max-alleles 2 $FILT_DIR/SNPs_indels_DP"$MIN_DP"_GQ"$MIN_GQ".tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/SNPs_indels_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz

# 3. Split SNPs and indels
bcftools filter -i "INFO/INDEL=1" $FILT_DIR/SNPs_indels_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz -O z --threads $CPU > $FILT_DIR/indels/indels_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
tabix -p vcf $FILT_DIR/indels/indels_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz

bcftools filter -e "INFO/INDEL=1" $FILT_DIR/SNPs_indels_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz -O z --threads $CPU > $FILT_DIR/SNPs/SNPs_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
#bcftools filter -i "INFO/INDEL=1" $DP_MISS_VCF -O z --threads $CPU > "$FILT_DIR/SNPs__"$MAX_ALL"all_maf"$MIN_MAF"_FM"$MAX_MISS"_minDP"$MIN_DP".vcf.gz"
tabix -p vcf $FILT_DIR/SNPs/SNPs_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
