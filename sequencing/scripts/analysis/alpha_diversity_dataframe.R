# look at alpha diversity 

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vegan)
library(gridExtra)
library(lme4)

# figure path
path_fig <- 'sequencing/plots/fresh'

# load data - latest run which we are happy with ####
# these files need to be there
ps <- readRDS('sequencing/data/output/20171024_17:18/ps_no_NA_phyla.rds')

# replace metadata with new metadata
# when wanting to add columns to metadata, it is better to edit metadata_creation and overwrite the metadata file as then it can be overwritten in all future files
meta_new <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE)
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 28,000. Woof.

# alpha diversity estimates ####

# prune OTUs that are not present in any of the samples
ps_sub <- prune_taxa(taxa_sums(ps) > 0, ps)

# metadata ###
m <- sample_data(ps_sub) %>%
  select(., SampleID, sample_name, treatment, evolution, preadapt_pop) %>%
  data.frame() %>%
  janitor::clean_names()

# calculate diversity measures of each sample ####
a_div <- estimate_richness(ps) %>%
  mutate(., SampleID = row.names(.)) %>%
  mutate(., pielou = Shannon / log(Observed)) %>%
  janitor::clean_names() %>%
  merge(., m, by = 'sample_id') %>%
  mutate_at(., c('sample_id', 'treatment', 'sample_name', 'evolution'), as.character) 

a_div <- select(a_div, sample_id, treatment, evolution, preadapt_pop, observed, pielou, inv_simpson) %>%
  filter(treatment != 'wt_ancestor')

# first plot of observed OTUs and evenness ####
gather(a_div, 'metric', 'value', c(observed, pielou, inv_simpson)) %>%
  ggplot(., aes(treatment, value)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  ggtitle('Diversity and evenness across treatments') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~metric, scale = 'free_y')

ggsave(file.path(path_fig, 'alpha_diversity.png'), last_plot(), height = 5, width = 10)

# look at individual clones only, does evolution with and without community impact diversity and evenness ####

# filter
a_div_ind <- filter(a_div, treatment == 'individual_clone')

# plot
ggplot(a_div_ind, aes(evolution, observed)) +
  MicrobioUoE::geom_pretty_boxplot(aes(fill = evolution, col = evolution)) +
  geom_point(aes(col = evolution), shape = 21, fill = 'white', position = position_jitter(width = 0.15)) +
  ggtitle('Observed OTUs') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# model
mod1 <- lmer(observed ~ evolution + (1|preadapt_pop), a_div_ind)
mod2 <- lmer(observed ~ 1 + (1|preadapt_pop), a_div_ind)
anova(mod1, mod2)

# pielous evenness
ggplot(a_div_ind, aes(evolution, pielou)) +
  MicrobioUoE::geom_pretty_boxplot(aes(fill = evolution, col = evolution)) +
  geom_point(aes(col = evolution), shape = 21, fill = 'white', position = position_jitter(width = 0.15)) +
  ggtitle("Pielou's evenness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mod1 <- lmer(pielou ~ evolution + (1|preadapt_pop), a_div_ind)
mod2 <- lmer(pielou ~ 1 + (1|preadapt_pop), a_div_ind)
anova(mod1, mod2)

# look at high diversity - 24 clones

# filter
a_div_24 <- filter(a_div, treatment %in% c('evolved_with_community', 'evolved_without_community'))

# plot
ggplot(a_div_24, aes(evolution, observed)) +
  MicrobioUoE::geom_pretty_boxplot(aes(fill = evolution, col = evolution)) +
  geom_point(aes(col = evolution), shape = 21, fill = 'white', position = position_jitter(width = 0.15)) +
  ggtitle('Observed OTUs') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# model
mod1 <- lm(observed ~ evolution, a_div_24)
mod2 <- lm(observed ~ 1, a_div_24)
anova(mod1, mod2)

# pielous evenness
ggplot(a_div_24, aes(evolution, pielou)) +
  MicrobioUoE::geom_pretty_boxplot(aes(fill = evolution, col = evolution)) +
  geom_point(aes(col = evolution), shape = 21, fill = 'white', position = position_jitter(width = 0.15)) +
  ggtitle("Pielou's evenness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mod1 <- lm(pielou ~ evolution, a_div_24)
mod2 <- lm(pielou ~ 1, a_div_24)
anova(mod1, mod2)

# look at diversity of pseudomonads across samples and treatments ####

SBW25 = "ACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAA"

ps2 <- subset_taxa(ps, Genus == 'Pseudomonas')
#ps2 <- subset_taxa(ps, rownames(tax_table(ps)) %in% c(SBW25))
sample_sums(ps2)

a_div_pseu <- #subset_taxa(ps2, rownames(tax_table(ps)) %in% c(SBW25)) %>%
  estimate_richness(ps2, measures = c('Shannon', 'Observed')) %>%
  mutate(., SampleID = row.names(.)) %>%
  mutate(., pielou = Shannon / log(Observed)) %>%
  janitor::clean_names() %>%
  merge(., m, by = 'sample_id') %>%
  mutate_at(., c('sample_id', 'treatment', 'sample_name', 'evolution'), as.character) 

a_div_pseu <- select(a_div_pseu, sample_id, treatment, evolution, preadapt_pop, observed, pielou) %>%
  filter(treatment != 'wt_ancestor')

gather(a_div_pseu, 'metric', 'value', c(observed, pielou)) %>%
  ggplot(., aes(treatment, value)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15, height = 0), size = 3) +
  ggtitle('Diversity and evenness across treatments') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~metric, scale = 'free_y')

# calculate max proportion of ASVs per sample
d_pseu <- psmelt(ps2) %>%
  janitor::clean_names() %>%
  filter(treatment != 'wt_ancestor') %>%
  mutate(SBW25 = ifelse(otu == SBW25, 'yes', 'no')) %>%
  group_by(sample, SBW25, treatment) %>%
  summarise(., abundance = sum(abundance)) %>%
  ungroup() %>%
  spread(., SBW25, abundance) %>%
  mutate(., prop = yes / (no + yes))

select(d_pseu, otu) %>%
  do(tibble(otu = paste('pseudomonad', 1:length(unique(.$otu)), sep = '_'), seq = unique(.$otu))) %>%
  write.csv(., 'sequencing/data/output/pseudomonads.csv', row.names = FALSE)
  
# plot
group_by(d_pseu, sample) %>%
  ggplot(., aes(treatment, prop)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  ggtitle('Proportion of SBW25 relative to the other Pseudomonads across treatments') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('proportion of reads assigned to SBW25') 

ggsave(file.path(path_fig, 'SBW25_success.pdf'), last_plot(), height = 5, width = 7)
ggsave(file.path(path_fig, 'SBW25_prop.png'), last_plot(), height = 5, width = 7)


# per OTU presence in samples
d_pres <- filter(d_pseu, abundance > 0) %>%
  group_by(otu) %>%
  summarise(., n_samples = n(),
            ave_prop = mean(prop))
