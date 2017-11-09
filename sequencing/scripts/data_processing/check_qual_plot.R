# raw processing ####

# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# Need to create the quality plots to examine the reads of each sample and set trimming parameters

# clean workspace before starting ####
rm(list = ls())

# set seed ####
set.seed(42)

# load in packages
library(dada2)
library(ggplot2)

# load in extra functions
source('scripts/Extra_Functions.R')

# list files in raw path ####
raw_fastq <- sort(list.files('data/raw_fastq', pattern = 'fast', full.names = TRUE, recursive = T))

# sort files for forward and reverse sequences ####
fnFs <- raw_fastq[grepl("R1", raw_fastq)]
fnRs <- raw_fastq[grepl("R2", raw_fastq)]

# Method 1 - create master Fwd and Rev Quality profiles ####
# create one big fastq file for Fwd and Rev and create single quality plots - suppose this could take a very long time for very big datasets

system(paste('cat', paste(fnFs, collapse = ' '), '> data/fwd_master.fastq', sep = ' '))
system(paste('cat', paste(fnFs, collapse = ' '), '> data/rev_master.fastq', sep = ' '))

# save plots out
pdf(file.path('plots', 'master_quality_profile.pdf'))
plotQualityProfile('data/fwd_master.fastq', n = 2e6) +
  ggtitle('Fwd reads master quality profile')
plotQualityProfile('data/rev_master.fastq', n = 2e6) +
  ggtitle('Rev reads master quality profile')
dev.off()

# Method 2 - single quality profile for every single sample ####
plot_qual(file.path('plots', 'qual_plot_raw_fastq.pdf'), fnFs, fnRs, height = 5, width = 7)

# examine this plot and set trimming parameters




