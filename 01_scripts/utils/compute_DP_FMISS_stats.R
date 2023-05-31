# Compute mean, min, max and sd of per-site total DP and F_MISSING, for raw and filtered data


# 1. Access file and import -----------------------------------------------
RAW_STATS_SNP <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/SNPs_DP_FMISS_tag.table'
FILT_STATS_SNP <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/SNPs/SNPs_MAF0.05_FMISS0.5_DP_FMISS_tag.table'

RAW_STATS_INDEL <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/indels_DP_FMISS_tag.table'
FILT_STATS_INDEL <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/indels/indels_MAF0.05_FMISS0.5_DP_FMISS_tag.table'

library(data.table)

raw_stats_SNP <- fread(RAW_STATS_SNP, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))
filt_stats_SNP <- fread(FILT_STATS_SNP, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))

raw_stats_indel <- fread(RAW_STATS_INDEL, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))
filt_stats_indel <- fread(FILT_STATS_INDEL, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))


# 2. Compute basic stats --------------------------------------------------
variants_stats <- function(x){
  # Compute mean
  print(paste('mean :', mean(x)))
  # Compute min and max
  print(paste('min :', min(x)))
  print(paste('max :', max(x)))
  # Compute sd
  print(paste('sd :', sd(x)))
}

# On raw, unfiltered data
variants_stats(raw_stats_SNP$DP)
variants_stats(raw_stats_SNP$F_MISSING)

variants_stats(raw_stats_indel$DP)
variants_stats(raw_stats_indel$F_MISSING)

# On filtered data
variants_stats(filt_stats_SNP$DP)
variants_stats(filt_stats_SNP$F_MISSING)

variants_stats(filt_stats_indel$DP)
variants_stats(filt_stats_indel$F_MISSING)
