#!/bin/bash


# manitou
# srun -c 2 -p small --time=1-00:00:00 --mem=50G -J indels_compare_FMISS_filters -o log/indels_compare_FMISS_filters_%j.log /bin/sh ./01_scripts/utils/compare_FMISS_filters_indels.sh &

# VARIABLES
# First 2 files were produced by the "utils/extract_unfiltered.sh" script
#RAW_SNP="/project/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/SNPs.vcf.gz"
RAW_INDEL="/project/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/indels.vcf.gz"

FILT_DIR="07_filtered"

REGIONS_EX="02_infos/excl_chrs.txt"

MIN_GQ=5
MAX_GQ=256 # deduct max allowed GQ from VCF, e.g. bcftools query -f "[%GQ\n]" $CALLS_DIR/raw/"$SAMPLE".vcf | sort -n | uniq | tail -n1

MIN_DP=4
MAX_DP=80 # arbitrary threshold of 5 * anticipated SR sequencing coverage, 5 * 16 

MIN_MAF=0.05
MAX_MAF=0.95

#MAX_MISS=0.5

CPU=2 

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load htslib/1.13


echo "Raw indels : $(zless $RAW_INDEL | grep -v ^'#' | wc -l) variants"

# 0. Remove non biallelic sites, assign missing genotype where DP and GQ are too low or too extreme, add tags and remove unwanted contigs from header
#bcftools view -m2 -M2 --threads $CPU $RAW_INDEL | bcftools +setGT --threads $CPU -- -t q -n . -e "FORMAT/DP >= $MIN_DP & FORMAT/DP < $MAX_DP & FORMAT/GQ >= $MIN_GQ & FORMAT/GQ < $MAX_GQ" | bcftools +fill-tags --threads $CPU -- -t AC,AC_Hom,AC_Het,AC_Hemi,AF,AN,NS,ExcHet,HWE,MAF,F_MISSING | bcftools sort | grep -vFf $REGIONS_EX > $FILT_DIR/indels/"$(basename -s .vcf.gz $RAW_INDEL)"_DP"$MIN_DP"_GQ"$MIN_GQ"_tagged.vcf
#bgzip $FILT_DIR/indels/"$(basename -s .vcf.gz $RAW_INDEL)"_DP"$MIN_DP"_GQ"$MIN_GQ"_tagged.vcf -f

#echo "With -m2 -M2, MIN_DP=$MIN_DP and MIN_GQ=$MIN_GQ : $(zless $FILT_DIR/indels/"$(basename -s .vcf.gz $RAW_INDEL)"_DP"$MIN_DP"_GQ"$MIN_GQ"_tagged.vcf.gz | grep -v ^'#' | wc -l) variants"


# Start with prefiltered SNPs and indels outputted by 03_filter_samples.sh
bcftools filter --threads $CPU -i "INFO/INDEL=1" $FILT_DIR/SNPs_indels_DP"$MIN_DP"_GQ"$MIN_GQ".vcf.gz | bcftools +fill-tags --threads $CPU -- -t AC,AC_Hom,AC_Het,AC_Hemi,AF,AN,NS,ExcHet,HWE,MAF,F_MISSING | bcftools sort -Oz > $FILT_DIR/indels/indels_DP"$MIN_DP"_GQ"$MIN_GQ"_tagged.vcf.gz

echo "With -m2 -M2, MIN_DP=$MIN_DP and MIN_GQ=$MIN_GQ : $(zless $FILT_DIR/indels/indels_DP"$MIN_DP"_GQ"$MIN_GQ"_tagged.vcf.gz | grep -v ^'#' | wc -l) indels"

# 1. Filter on MAF and proportion of missing genotypes, with F_MISS = 0.75, 0.625, 0.125, 0.25, 0.375 

for MAX_MISS in 0.75 0.625 0.5 0.375 0.25 0.125 ; 
do 
  #bcftools view --max-alleles 2 $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_DP"$MIN_DP"_GQ"$MIN_GQ"_tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_DP"$MIN_DP"_GQ"$MIN_GQ"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
  bcftools view --max-alleles 2 $FILT_DIR/indels/indels_DP"$MIN_DP"_GQ"$MIN_GQ"_tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/indels/indels_DP"$MIN_DP"_GQ"$MIN_GQ"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
  echo "with $MAX_MISS : $(zless $FILT_DIR/indels/indels_DP"$MIN_DP"_GQ"$MIN_GQ"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants" ;
done


#MAX_MISS=0.75
#bcftools view --max-alleles 2 $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
#echo "with $MAX_MISS : $(zless $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"


#MAX_MISS=0.625
#bcftools view --max-alleles 2 $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
#echo "with $MAX_MISS : $(zless $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"

#MAX_MISS=0.5
#bcftools view --max-alleles 2 $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
#echo "with $MAX_MISS : $(zless $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"

#MAX_MISS=0.375
#bcftools view --max-alleles 2 $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
#echo "with $MAX_MISS : $(zless $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"


#MAX_MISS=0.25
#bcftools view --max-alleles 2 $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
#echo "with $MAX_MISS : $(zless $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"


#MAX_MISS=0.125
#bcftools view --max-alleles 2 $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
#echo "with $MAX_MISS : $(zless $FILT_DIR/"$(basename -s .vcf.gz $RAW_INDEL)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants"

