# Compute mean, min, max and sd of per-site total DP and F_MISSING, for raw and filtered data


# 1. Access file and import -----------------------------------------------
RAW_STATS_SNP <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/SNPs_DP_FMISS_tag.table'
FILT_STATS_SNP <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/SNPs/SNPs_MAF0.05_FMISS0.5_DP_FMISS_tag.table'

RAW_STATS_INDEL <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/indels_DP_FMISS_tag.table'
FILT_STATS_INDEL <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/indels/indels_MAF0.05_FMISS0.5_DP_FMISS_tag.table'

library(data.table)
library(dplyr)
library(ggplot2)


raw_stats_SNP <- fread(RAW_STATS_SNP, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))
filt_stats_SNP <- fread(FILT_STATS_SNP, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))

raw_stats_indel <- fread(RAW_STATS_INDEL, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))
filt_stats_indel <- fread(FILT_STATS_INDEL, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))


# 2. Compute basic stats --------------------------------------------------
variants_stats <- function(x){
  # Total 
  print(paste('number of variants :', length(x)))
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

# 3. Explore F_MISSING ----------------------------------------------------
# Split F_MISSING into equal size bins
raw_stats_SNP$F_MISS_bins <- cut_interval(raw_stats_SNP$F_MISSING, length = 0.1, right = FALSE)

# Assign 'kept' or 'filtered out' tag to each SNP
raw_stats_SNP$F_MISS_groups <-
  sapply(X = raw_stats_SNP$F_MISS_bins,
         FUN = function(x){
           ifelse(x %in% levels(raw_stats_SNP$F_MISS_bins)[1:5], 
                  yes = 'passed',
                  no = 'failed')
         })


# Plot
## add type in case we need it for combined plotting
raw_stats_SNP$TYPE <- 'SNPs'
write.table(raw_stats_SNP,
            col.names = c('CHROM', 'POS', 'DP', 'F_MISSING',
                          'F_MISS_bins', 'F_MISS_groups', 'TYPE'),
            file = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/SNPs/SNPs_MAF0.05_FMISS0.5_DP_FMISS_FMISS_stats.table',
            sep = '\t', quote = FALSE, row.names = FALSE)

SNPs_F_MISS_plot <-
ggplot(data = raw_stats_SNP) +
  geom_bar(aes(F_MISS_bins, fill = F_MISS_groups)) + 
  theme(
    axis.text.x = element_text(angle = 45, size = 8, hjust = 1)
  ) + scale_fill_manual(values = c('red', 'grey60')) +
  labs(x = 'Proportion of missing genotypes',
       y = 'Count',
       fill = 'F_MISSING filter')

saveRDS(SNPs_F_MISS_plot,
        file = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/SNPs/SNPs_MAF0.05_FMISS0.5_DP_FMISS_FMISSplot.rds')



# Split F_MISSING into equal size bins
raw_stats_indel$F_MISS_bins <- cut_interval(raw_stats_indel$F_MISSING, length = 0.1, right = FALSE)

# Assign 'kept' or 'filtered out' tag to each SNP
raw_stats_indel$F_MISS_groups <-
  sapply(X = raw_stats_indel$F_MISS_bins,
         FUN = function(x){
           ifelse(x %in% levels(raw_stats_indel$F_MISS_bins)[1:5], 
                  yes = 'passed',
                  no = 'failed')
         })


# Plot
## add type in case we need it for combined plotting
raw_stats_indel$TYPE <- 'Indels'
write.table(raw_stats_indel,
            col.names = c('CHROM', 'POS', 'DP', 'F_MISSING',
                          'F_MISS_bins', 'F_MISS_groups', 'TYPE'),
            file = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/indels/indels_MAF0.05_FMISS0.5_DP_FMISS_FMISS_stats.table',
            sep = '\t', quote = FALSE, row.names = FALSE)



indels_F_MISS_plot <-
ggplot(data = raw_stats_indel) +
  geom_bar(aes(F_MISS_bins, fill = F_MISS_groups)) + 
  theme(
    axis.text.x = element_text(angle = 45, size = 8, hjust = 1)
  ) + scale_fill_manual(values = c('red', 'grey60')) +
  labs(x = 'Proportion of missing genotypes',
       y = 'Count',
       fill = 'F_MISSING filter')

saveRDS(indels_F_MISS_plot,
        file = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/indels/indels_MAF0.05_FMISS0.5_DP_FMISS_FMISSplot.rds')



