# first actual analyses ####

# clear workspace ####
rm(list = ls())

# load packages ####
library(phyloseq)
library(tidyverse)
library(magrittr)

# functions to use ####
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}

# load in data ####
ps <- readRDS('DNA/data/output/20170125_16:43/20170125_16:43_ps.rds')

# alpha diversity estimates ####
# prune OTUs that are not present in any of the samples
ps_sub <- prune_taxa(taxa_sums(ps) > 0, ps)

# prune for just autotrophs and cyanobacteria
ps_sub <- subset_taxa(ps, Kingdom == 'Eukaryota' | Phylum == 'Cyanobacteria')

# first plot
plot_richness(ps_sub, x = 'ancestral', col = 'treatment', measures = 'Chao1') +
  geom_boxplot() +
  scale_colour_manual(values = c('blue', 'green4', 'red'))

# calculate diversity measures of each sample ####
a_diversity <- estimate_richness(ps_sub) %>%
  mutate(., SampleID = mgsub(c('X', '\\.'), c('', '-'), row.names(.)))

# metadata ###
m <- sample_data(ps_sub) %>%
  select(., SampleID, sample, treatment, ancestral) %>%
  data.frame()

a_diversity <- merge(a_diversity, m, by = 'SampleID') %>%
  mutate_at(., c('sample', 'treatment', 'ancestral'), as.character) %>%
  mutate(time_point = ifelse(treatment == 'T0', '0', '1'),
         pond = ifelse(grepl('C', sample), 'C', sample),
         ancestral = case_when(.$ancestral == 'comb' & .$treatment == 'A' ~ 'comb_A',
                               .$ancestral == 'comb' & .$treatment == 'W' ~ 'comb_W',
                               .$sample == 'CA' ~ 'comb_A',
                               .$sample == 'CW' ~ 'comb_W',
                               TRUE ~ .$ancestral))

# look at T0 differences
ggplot(filter(a_diversity)) +
  geom_boxplot(aes(ancestral, Observed, col = treatment)) +
  scale_colour_manual(values = c('blue', 'black', 'red')) +
  facet_wrap(~ ancestral, scales = 'free_x')


# duplicate T0 rows for point plotting
# duplicate pond T0 samples once 
T0 <- rbind(filter(a_diversity, treatment == 'T0' & !sample %in% c('CW', 'CA')), filter(a_diversity, treatment == 'T0' & !sample %in% c('CW', 'CA'))) %>%
  mutate(., treatment = ifelse(row_number() <= n()/2, 'W', 'A'))

# Duplicate combined samples ten times
C_T0 <- filter(a_diversity, treatment == 'T0' & sample %in% c('CW', 'CA'))
C_T0 <- C_T0[rep(row.names(C_T0), each = 10),] %>%
  mutate(., n = rep(1:10, times = 2),
         treatment = ifelse(ancestral == 'comb_W', 'W', 'A')) %>%
  select(., -sample) %>%
  unite(., sample, c(pond, n), sep = '', remove = FALSE) %>%
  select(., -n)

# bind together  
a_div <- bind_rows(filter(a_diversity, treatment != 'T0'), T0, C_T0) %>%
  mutate(., treatment = plyr::mapvalues(treatment, from = c("A", 'W'), to = c('amb', 'warm')),
         AncIsTreat = ifelse(treatment == ancestral, 'Yes', 'No'))

# plot lines ###
ggplot(a_div) +
  geom_point(aes(time_point, Shannon, col = AncIsTreat)) +
  geom_line(aes(time_point, Shannon, group = sample)) +
  scale_colour_manual(values = c('blue', 'red')) +
  facet_wrap(~ treatment)

ggplot(a_div) +
  geom_point(aes(time_point, Shannon, col = AncIsTreat)) +
  geom_line(aes(time_point, Shannon, col = AncIsTreat, group = treatment)) +
  scale_colour_manual(values = c('dark grey', 'black')) +
  facet_wrap(~ sample + ancestral, labeller = labeller(.multi_line = FALSE)) +
  theme_bw()

# plot observed at T0 and TFinal
ggplot(a_div) +
  geom_point(aes(time_point, Observed, col = ancestral)) +
  geom_line(aes(time_point, Observed, col = ancestral, group = sample)) +
  scale_colour_manual(values = c('blue', 'black', 'black', 'red')) +
  facet_wrap(~ treatment, labeller = labeller(.multi_line = FALSE)) +
  theme_bw()

select(a_div, Observed, sample, treatment, ancestral, pond, AncIsTreat, time_point) %>%
  spread(., time_point, Observed) %>%
  rename(., start = `0`, end = `1`) %>%
  ggplot(.) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(start, end, col = treatment), position = position_jitter(width = 0.5)) +
  theme_bw() +
  xlab('Number of OTUs at start') +
  ylab('Number of OTUs at end') +
  scale_colour_manual(values = c('blue', 'red'))

