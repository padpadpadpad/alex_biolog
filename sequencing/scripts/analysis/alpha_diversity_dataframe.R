# look at alpha diversity 

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vegan)
library(gridExtra)

# function for creating the box plots I like
add_boxplot <- function(...){
  list(
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), col = 'black', ...),
    stat_summary(geom = 'crossbar', position = position_dodge(width = 0.75), fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))})
  )
}

# figure path
path_fig <- 'sequencing/plots'

# load data - latest run which we are happy with ####
ps <- readRDS('sequencing/data/output/20171024_17:18/20171024_17:18_ps.rds')

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 30,000. Woof.

# alpha diversity estimates ####

# prune OTUs that are not present in any of the samples
ps_sub <- prune_taxa(taxa_sums(ps) > 0, ps)

# metadata ###
m <- sample_data(ps_sub) %>%
  select(., SampleID, sample_name, treatment) %>%
  data.frame()

# calculate diversity measures of each sample ####
a_div <- estimate_richness(ps) %>%
  mutate(., SampleID = row.names(.)) %>%
  merge(., m, by = 'SampleID') %>%
  mutate_at(., c('SampleID', 'treatment', 'sample_name'), as.character) 

# a few plots ####
ggplot(a_div, aes(treatment, Observed)) +
  add_boxplot(fill = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15)) +
  ggtitle('Observed OTUs')
ggplot(a_div, aes(treatment, Shannon)) +
  add_boxplot(fill = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15)) +
  ggtitle('Shannon diversity')
ggplot(a_div, aes(treatment, Fisher)) +
  add_boxplot(fill = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15)) +
  ggtitle('Fishers diversity')
ggplot(a_div, aes(treatment, Chao1)) +
  add_boxplot(fill = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15)) +
  ggtitle('Chao index')


