#!/bin/bash

# Filter SNPs 
# Run this script from SNPs_indels_SR directory using this command :
# srun -c 4 --mem=10G -p ibis_medium --time=7-00:00:00 -J 03_filter -o log/03_filter_%j.log 01_scripts/03_filter.sh &

# VARIABLES
CHR_LIST="02_infos/chr.txt"
BAM_LIST="02_infos/bam_list.txt"

CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
MERGED_VCF="$MERGED_DIR/SNPs_indels.vcf.gz"
VCF_LIST="02_infos/VCF_list.txt"

CPU=4

## Filtering criteria : 
MIN_DP=1      # min depth per site (number of reads overlapping given site)
MAX_ALL=2     # max number of alleles per site
MIN_MAF=0.05  # min alt allele frequency
MAX_MAF=0.95  # max alt allele frequency
MAX_MISS=0.5  # max proportion of missing genotypes allowed (fraction of total sample count)

TAGGED_VCF="$MERGED_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_tagged.vcf.gz"
MAX_ALL_VCF="$FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_"$MAX_ALL"all.vcf.gz"
MIN_MAF_VCF="$FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_"$MAX_ALL"all_maf"$MIN_MAF".vcf.gz"
DP_MISS_VCF="$FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_"$MAX_ALL"all_maf"$MIN_MAF"_FM"$MAX_MISS"_minDP"$MIN_DP".vcf.gz"

# LOAD REQUIRED MODULES
module load bcftools

# 1. Add tags
bcftools +fill-tags $MERGED_VCF -Oz -- -t all --threads $CPU > $TAGGED_VCF

# Filter with same criteria as SVs
# 2. Filter for max number of alleles = extract biallelic sites
bcftools view --max-alleles $MAX_ALL $TAGGED_VCF -O z --threads $CPU > $MAX_ALL_VCF

# 3. Filter for maf  
bcftools filter -i "INFO/MAF >= $MIN_MAF" $MAX_ALL_VCF -O z --threads $CPU > $MIN_MAF_VCF

# 4. Filter for depth, where at least 50% of samples have been genotyped
#bcftools filter --threads $CPU -S . -e "FORMAT/DP > $MIN_DP | FORMAT/DP < $MAX_DP" $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_"$MAX_ALL"all_maf"$MIN_MAF".vcf.gz -Ou | 
#  bcftools view --threads $CPU -i "F_MISSING < $MAX_MISS " -Oz -o $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_"$MAX_ALL"all_maf"$MIN_MAF"_FM"$MAX_MISS"_minDP"$MIN_DP".vcf.gz
#bcftools filter -i "N_PASS(GT != 'mis' & FMT/DP > $MIN_DP) > 30"

# Set genotypes to missing if sample's DP < $MIN_DP, then filter for sites genotyped in > $MAX_MISS proportion of samples
bcftools +setGT $MIN_MAF_VCF -- -t q -n . -i "FORMAT/DP < $MIN_DP" --threads $CPU | bcftools view -i "F_MISSING < $MAX_MISS " -O z --threads $CPU > $DP_MISS_VCF

# 5. Index with
tabix -p vcf $DP_MISS_VCF

# 6. Split SNPs and indels
bcftools filter -i "INFO/INDEL=1" $DP_MISS_VCF -O z --threads $CPU > $FILT_DIR/filtered_indels.vcf.gz
#bcftools filter -i "INFO/INDEL=1" $DP_MISS_VCF -O z --threads $CPU > "$FILT_DIR/indels__"$MAX_ALL"all_maf"$MIN_MAF"_FM"$MAX_MISS"_minDP"$MIN_DP".vcf.gz"

tabix -p vcf $FILT_DIR/filtered_indels.vcf.gz

bcftools filter -e "INFO/INDEL=1" $DP_MISS_VCF -O z --threads $CPU > $FILT_DIR/filtered_SNPs.vcf.gz
#bcftools filter -i "INFO/INDEL=1" $DP_MISS_VCF -O z --threads $CPU > "$FILT_DIR/SNPs__"$MAX_ALL"all_maf"$MIN_MAF"_FM"$MAX_MISS"_minDP"$MIN_DP".vcf.gz"

tabix -p vcf $FILT_DIR/filtered_SNPs.vcf.gz
