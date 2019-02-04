# analysis of whether pairwise distance from biolog plate correlates with distance from effect on community composition

# analysis of number of clone in sample ####
rm(list = ls())

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(viridis)
library(tibble)
library(patchwork)
library(widyr)

#--------------#
# functions ####
#--------------#

# convert distance matrix to dataframe
dist_2_df <- function(dist_ob){
  m <- as.matrix(dist_ob) # coerce dist object to a matrix
  xy <- t(combn(colnames(m), 2))
  return(data.frame(xy, dist=m[xy], stringsAsFactors = FALSE))
}

#-------------------------------------#
# setup workspace ####
#-------------------------------------#

# set seed
set.seed(42)

# figure path
path_fig <- 'plots'

#--------------------------------#
# prepare sequencing analysis ####
#--------------------------------#

# load data - latest run which we are happy with
ps <- readRDS('sequencing/data/output/20171024_17:18/ps_no_NA_phyla.rds')

# replace metadata with new metadata
meta_new <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE)
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 28,000 Woof.

# initial subsetting - do not want nmc_t0 or wt_ancestor
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0', 'wt_ancestor'))
ps2 <- prune_samples(to_keep$SampleID, ps)

# remove Pseudomonas fluorescens reads from the data, leave other pseudomonads in
SBW25 = "ACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAA"

ps2 <- subset_taxa(ps2, ! rownames(tax_table(ps2)) %in% c(SBW25))

# samples to keep - just individual clones
to_keep <- filter(meta_new, treatment %in% c('individual_clone'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
# weighted Unifrac - relative abundances
# Bray-Curtis - absolute abundances
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution),
                 preadapt_pop = as.factor(preadapt_pop),
                 pop2 = preadapt_pop) %>%
  column_to_rownames(., 'SampleID')

# convert pre_adapt pop to letters
levels(d_samp$pop2) <- letters[1:length(levels(d_samp$pop2))]

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac', upper = TRUE)

# replace these labels with the sample_name
attr(ps_wunifrac, 'Labels') <- d_samp$sample_name

# convert to dataframe 
d_wunifrac <- dist_2_df(ps_wunifrac)

#------------------------#
# prepare biolog data ####
#------------------------#

# read in all files
d <- read.csv('biolog/data/20181203_processed.csv', stringsAsFactors = FALSE)

# filter out blank
blank <- filter(d, sample == 'M9')

blank_ave <- group_by(blank, substrate, od_wave) %>%
  summarise(., blank = median(od),
            sd_blank = sd(od))

d <- filter(d, sample != 'M9') %>%
  merge(., blank_ave, by = c('substrate', 'od_wave')) %>%
  mutate(., od_cor = od - blank)

# lets use T4 590
d_t4_590 <- filter(d, od_wave == 590 & tp == 'T4')

# remove obvious outliers
d_t4_590 <- filter(d_t4_590, sample != 'anc1')
d_t4_590 <- filter(d_t4_590,! sample %in% c('33', '45'))

# add column for ranked mean OD_cor
d_t4_590 <- group_by(d_t4_590, substrate) %>%
  mutate(mean_od = mean(od_cor)) %>%
  ungroup() %>%
  mutate(., rank = dense_rank(desc(mean_od)))

# filter so that only wells where growth could occur were considered
d_t4_590 <- filter(d_t4_590, mean_od > 0.05)

# find pairwise cor between clones 
pop_dists_df <- filter(d_t4_590, evolved %in% c('with_community', 'without_community')) %>%
  pairwise_dist(sample, substrate, od_cor, upper = TRUE) %>%
  rename(clone_i = item1, clone_j = item2)

head(pop_dists_df)

# combine tables euclidean distance and weighted Unifrac distance matrices together
d <- merge(rename(d_wunifrac, clone_i = X1, clone_j = X2), pop_dists_df, , by = c('clone_i', 'clone_j'))

d <- rename(d, euclid_dist = distance, wunifrac_dist = dist) %>%
  merge(., select(d_samp, clone_i = sample_name, evol_i = evolution, pop_i = preadapt_pop), by = 'clone_i') %>%
  merge(., select(d_samp, clone_j = sample_name, evol_j = evolution, pop_j = preadapt_pop), by = 'clone_j') %>%
  mutate(., same_pop = ifelse(pop_i == pop_j, 'sympatry', 'allopatry'),
         same_treat = ifelse(evol_i == evol_j, 'same_treat', 'oppo_treat'))

ggplot(d, aes(euclid_dist, wunifrac_dist, col = same_treat)) +
  geom_point() + 
  facet_wrap(~ evol_i)

ggplot(d, aes(euclid_dist, wunifrac_dist, col = evol_i)) +
  geom_point()

# no correlation, but something else
d2 <- merge(d_wunifrac, select(d_samp, X1 = sample_name, evol_1 = evolution, pop_1 = preadapt_pop), by = 'X1') %>%
  merge(., select(d_samp, X2 = sample_name, evol_2 = evolution, pop_2 = preadapt_pop), by = 'X2') %>%
  mutate(., same_pop = ifelse(pop_1 == pop_2, 'sympatry', 'allopatry'),
         same_treat = ifelse(evol_1 == evol_2, 'same_treat', 'oppo_treat')) %>%
  gather(., 'random', 'clone', c(X1, X2)) %>%
  group_by(clone, same_treat) %>%
  summarise(., ave_dist = max(dist)) %>%
  ungroup()

head(d2)

d2 <- merge(d2, select(d_samp, clone = sample_name, evolution, preadapt_pop), by = 'clone')

ggplot(d2, aes(same_treat, ave_dist)) +
  geom_pretty_boxplot(fill = 'black', col = 'black')
