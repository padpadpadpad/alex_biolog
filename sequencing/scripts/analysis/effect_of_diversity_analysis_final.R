# analysis of number of clone in sample ####
rm(list = ls())

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vegan)
library(gridExtra)
library(viridis)
library(ggridges)
library(tibble)
library(patchwork) # devtools::install_github('thomasp85/patchwork')
library(lme4)
# if not installed, install mctoolsr run devtools::install_github('leffj/mctoolsr')

#--------------#
# functions ####
#--------------#

# code stolen from phyloseq website
get_top_taxa <- function(ps, tax_rank, to_keep){
  temp <- tapply(phyloseq::taxa_sums(ps), phyloseq::tax_table(ps)[, tax_rank], sum, na.rm = TRUE)
  temp2 <-  names(sort(temp, TRUE))[1:to_keep]
  return(temp2)
}

# get betadisper dataframes

# getting distances from betadisper() object
betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}

# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# betadisper data
get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}

# get distance from 00
dist_between_points <- function(x1, x2, y1, y2){
  return(abs(sqrt((x1 - x2)^2+(y1-y2)^2)))
}

#-------------------------------------#
# setup workspace and load in data ####
#-------------------------------------#

# set seed
set.seed(42)

# figure path
path_fig <- 'plots'

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

# -------------------------------------------------------------------------------#
# Analysis of the preadaptation with and without nmc on community composition ####
#--------------------------------------------------------------------------------#

# perfectly replicated for 1, 4 and 24 clones
# haver to do these separately because of nesting in the individual clone treatment

# 1. individual clones ####

# samples to keep - just individual clones
to_keep <- filter(meta_new, treatment %in% c('individual_clone'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution),
                 preadapt_pop = as.factor(preadapt_pop)) %>%
  column_to_rownames(., 'SampleID')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# get ready to permute distances between centroids
vec_1 <- tibble::rownames_to_column(d_samp, var = 'sample') %>%
  filter(., evolution == 'with_community') %>% 
  pull(sample)
vec_2 <- tibble::rownames_to_column(d_samp, var = 'sample') %>%
  filter(., evolution == 'without_community') %>% 
  pull(sample)

# run an Adonis test
# tried using strata to shuffle within but not across replicates within each treatment
# not right atm, but no effect anyway
mod_div_1 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_1 <- betadisper(ps_wunifrac, d_samp$evol_fac)
d_perm1 <- permute_distance(ps_wunifrac, vec_1, vec_2) %>%
  mutate(., nclones = 1)

# 2. 4 related clones ####
# filter treatments to keep
to_keep <- filter(meta_new, treatment %in% c('4_related_clones'))
ps_sub <- prune_samples(to_keep$SampleID, ps)

# remove Pseudomonas reads
ps_sub <- subset_taxa(ps_sub, Genus != 'Pseudomonas')

# transform counts to relative abundances for ordination
# cannot really find that much consenses on either side of the argument
# whether or not to do this appears to be a little controversial still
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution)) %>%
  column_to_rownames(., 'SampleID')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# get ready to permute distances between centroids
vec_1 <- tibble::rownames_to_column(d_samp, var = 'sample') %>%
  filter(., evolution == 'with_community') %>% 
  pull(sample)
vec_2 <- tibble::rownames_to_column(d_samp, var = 'sample') %>%
  filter(., evolution == 'without_community') %>% 
  pull(sample)

usedist::dist_between_centroids(ps_wunifrac, vec_1, vec_2)

# run an Adonis test
mod_div_4 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_4 <- betadisper(ps_wunifrac, d_samp$evol_fac)
d_perm4 <- permute_distance(ps_wunifrac, vec_1, vec_2) %>%
  mutate(., nclones = 4)

# 3. 24 related clones ####

# filter samples to keep
to_keep <- filter(meta_new, treatment %in% c('evolved_with_community', 'evolved_without_community'))
ps_sub <- prune_samples(to_keep$SampleID, ps)

# remove Pseudomonas reads
ps_sub <- subset_taxa(ps_sub, Genus != 'Pseudomonas')

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# do an adonis
# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution)) %>%
  column_to_rownames(., 'SampleID')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# get ready to permute distances between centroids
vec_1 <- tibble::rownames_to_column(d_samp, var = 'sample') %>%
  filter(., evolution == 'with_community') %>% 
  pull(sample)
vec_2 <- tibble::rownames_to_column(d_samp, var = 'sample') %>%
  filter(., evolution == 'without_community') %>% 
  pull(sample)

# run an Adonis test
mod_div_24 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_24 <- betadisper(ps_wunifrac, d_samp$evol_fac)
d_perm24 <- permute_distance(ps_wunifrac, vec_1, vec_2) %>%
  mutate(., nclones = 24)

# Figure 1, looking at effect of pre-adaptation context across different levels of diversity ####

# Make big plot and then grab the data
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0', '4_unrelated_clones',  'wt_ancestor'))

ps_sub <- prune_samples(to_keep$SampleID, ps)

# remove Pseudomonas reads
ps_sub <- subset_taxa(ps_sub, Genus != 'Pseudomonas')

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution)) %>%
  column_to_rownames(., 'SampleID') %>%
  unite(., 'id', nclones_fac, evol_fac, sep = ':')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# run a betadisper
mod_betadisper <- betadisper(ps_wunifrac, d_samp$id)

# grab centroids and other data
d_fig1 <- get_betadisper_data(mod_betadisper)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(d_fig1$centroids, group, PCoA1, PCoA2), select(d_fig1$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
d_fig1$eigenvector$distances <- d_fig1$distances$distances

# split up group into clones and evolution context
betadisper_lines <- separate(betadisper_lines, group, c('nclones', 'evol'), sep =':') %>%
  mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                             nclones == 'C_high' ~ 'negative_control',
                             TRUE ~ evol),
         nclones =forcats::fct_relevel(nclones,
                                   "C_1", "C_4", "C_24"))
d_fig1$centroids <- separate(d_fig1$centroids, group, c('nclones', 'evol'), sep =':') %>%
  mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                             nclones == 'C_high' ~ 'negative_control',
                             TRUE ~ evol),
         nclones =forcats::fct_relevel(nclones,
                                       "C_1", "C_4", "C_24"))
d_fig1$eigenvector <- separate(d_fig1$eigenvector, group, c('nclones', 'evol'), sep =':') %>%
  mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                             nclones == 'C_high' ~ 'negative_control',
                             TRUE ~ evol),
         nclones =forcats::fct_relevel(nclones,
                                       "C_1", "C_4", "C_24"))

# clone labels
labels = c('1 clone', '4 clones', '24 clones')
facet <- c(C_1 = '1 clone', C_4 = '4 clones', C_24 = '24 clones')

# plot PCoA
fig1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, alpha = 1 - distances), select(filter(d_fig1$eigenvector, evol %in% c('negative_control')), -nclones), size = 0.75, col = 'blue') +
  geom_point(aes(PCoA1, PCoA2), select(filter(d_fig1$centroids, evol %in% c('negative_control')), -nclones), size = 5, col = 'blue') +
  geom_point(aes(PCoA1, PCoA2, alpha = 1 - distances), select(filter(d_fig1$eigenvector, evol %in% c('lacz_ancestor')), -nclones), size = 0.75, col = 'orange') +
  geom_point(aes(PCoA1, PCoA2), select(filter(d_fig1$centroids, evol %in% c('lacz_ancestor')), -nclones), size = 5, col = 'orange') +
  geom_point(aes(PCoA1, PCoA2, col = evol, alpha = 1 - distances), filter(d_fig1$eigenvector, ! evol %in% c('lacz_ancestor', 'negative_control')), size = 0.75) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(filter(betadisper_lines, ! evol %in% c('lacz_ancestor', 'negative_control'))), col = evol, alpha = 1 - distances), filter(betadisper_lines, ! evol %in% c('lacz_ancestor', 'negative_control'))) +
  geom_point(aes(PCoA1, PCoA2, col = evol), filter(d_fig1$centroids, ! evol %in% c('lacz_ancestor', 'negative_control')), size = 5) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ylab('PCoA Axis 2') +
  xlab('PCoA Axis 1') +
  theme(legend.position = 'bottom') +
  facet_wrap(~ nclones, labeller = labeller(nclones = facet)) +
  ggtitle('PCoA of the effect of preadaptation history across different levels of diversity',
          subtitle = 'blue points are the nmc, orange are lacz ancestor') +
  scale_alpha(range = c(0.0001, 0.75), guide = FALSE) +
  scale_color_manual('', values = c('grey', 'black'))

# save plot, other ways are available
ggsave(file.path(path_fig, 'effect_of_evol_history.png'), fig1, height = 6, width = 12)

# plot all distances between centroids
d_perm <- bind_rows(d_perm1, d_perm4, d_perm24) %>%
  group_by(., nclones, type) %>%
  tidybayes::mean_qi() %>%
  ungroup() %>%
  mutate(., nclones_fac = paste('C_', nclones, sep = ''))

fig2 <- ggplot() +
  geom_point(aes(forcats::fct_reorder(nclones_fac, nclones), dist), filter(d_perm, type == 'actual_distance'), size = 3) +
  geom_linerange(aes(forcats::fct_reorder(nclones_fac, nclones), ymin = .lower, ymax = .upper), filter(d_perm, type == 'permuted')) +
  ylim(c(0, 0.15)) +
  theme_bw() +
  ylab('pairwise distance between centroids') +
  xlab('clonal diversity') +
  scale_x_discrete(labels = labels)

fig3 <- fig1 + fig2 + plot_layout(ncol = 2, widths = c(0.75, 0.25))

ggsave(file.path(path_fig, 'effect_of_evol_history.png'), fig3, height = 4.5, width = 12)

#-----------------------------------------------------------------------#
# Weighted Unifrac just on diversity, ignoring preadaptation context ####
#-----------------------------------------------------------------------#

# change some values of nclones
meta_new <- sample_data(ps) %>% data.frame() %>%
  mutate(., n_clones = paste('C_', n_clones, sep = ''))
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# filter samples
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0',  'wt_ancestor'))
ps_sub <- prune_samples(to_keep$SampleID, ps)

# remove Pseudomonas reads
ps_sub <- subset_taxa(ps_sub, Genus != 'Pseudomonas')

# make counts proportions
ps_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# get the distance matrix out of the data
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp,
                 nclones_fac = as.factor(n_clones)) %>%
  column_to_rownames(., 'SampleID')

# run an Adonis test
mod_nclonefac <-  vegan::adonis(ps_wunifrac ~ nclones_fac, data = d_samp, n_perm = 9999)

# run a multiple comparison to see which treatments are different
mult_comp <- mctoolsr::calc_pairwise_permanovas(ps_wunifrac, d_samp, 'nclones_fac', n_perm = 9999) %>%
  mutate(., pvalHolm = p.adjust(pval, method = 'holm'))

# save dataset out
saveRDS(mult_comp, 'sequencing/data/output/mult_comp.rds')
# loads of comparisons. Significant differences will be determined by p value correction.

# overwrite metadata to allow plotting of nclones_fac
sample_data(ps_prop) <- sample_data(d_samp)

# beta-diversity analysis - look at homogeneity of variances
mod1_dispers <- betadisper(ps_wunifrac, d_samp$nclones_fac)

# plot of model
plot(mod1_dispers)
boxplot(mod1_dispers)

# anova
anova(mod1_dispers)

# Permutation test for F
pmod <- permutest(mod1_dispers, pairwise = TRUE)

# Tukey's Honest Significant Differences
T_HSD <- TukeyHSD(mod1_dispers)

# plot distances to centroid - try and play with alpha

# get betadisper data ####
betadisper_dat <- get_betadisper_data(mod1_dispers)

# do some transformations on the data
betadisper_dat$eigenvalue <- mutate(betadisper_dat$eigenvalue, percent = eig/sum(eig)*100)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(betadisper_dat$centroids, group, PCoA1, PCoA2), select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
betadisper_dat$eigenvector$distances <- betadisper_dat$distances$distances

# plot PCoA
p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group, alpha = 1 - distances), betadisper_dat$eigenvector, size = 0.75) +
  #geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull, alpha = 0.1) +
  #geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = group, alpha = 1 - distances), betadisper_lines) +
  stat_ellipse(aes(PCoA1, PCoA2, col = group, group = group), type = "t", betadisper_dat$eigenvector) +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 5) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('PCoA Axis 2') +
  xlab('PCoA Axis 1') +
  theme(legend.position = 'bottom') +
  coord_fixed(sqrt(betadisper_dat$eigenvalue$percent[2]/betadisper_dat$eigenvalue$percent[1])) +
  ggtitle('PCoA plot looking at the effect of diversity') +
  scale_alpha(range = c(0.0001, 0.75), guide = FALSE)

# plot distances from centroid
p2 <- ggplot(betadisper_dat$distances, aes(forcats::fct_relevel(group, 'C_24', after = 2), distances, fill = group, col = group)) +
  geom_boxplot(aes(fill = group, col = group), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(group, distances, col = group), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.3)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('') +
  ggtitle('Distance to centroid across diversity')

p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(2.5, 1))

ggsave(file.path(path_fig, 'PCoA_plot_diversity.png'), p3, height = 5, width = 13)

# Nothing to be used below this!!!! ####

# filter out individual clones and see if fitness predicts distance from lacz ancestor

# filter samples
to_keep <- filter(meta_new, n_clones %in% c('1', 'high') & treatment != 'wt_ancestor')
ps_sub <- prune_samples(to_keep$SampleID, ps)

# remove Pseudomonas reads
ps_sub <- subset_taxa(ps_sub, Genus != 'Pseudomonas')

# make counts proportions
ps_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# ordinate and plot
ord_wUni <- ordinate(ps_prop, method = 'MDS', distance = 'wunifrac')

evals <- ord_wUni$values$Eigenvalues

# plot
plot_ordination(ps_prop, ord_wUni, color = "fitness") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~ treatment)

# get weighted-unifrac matrix
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp,
                 nclones_fac = as.factor(n_clones)) %>%
  column_to_rownames(., 'SampleID')

d_wunifrac <- as.matrix(ps_wunifrac, labels = TRUE) %>%
  as.data.frame.table(., stringsAsFactors = FALSE)

d_meta_ind_clone <- sample_data(ps_prop) %>% 
  data.frame() %>%
  filter(., treatment == 'individual_clone') %>%
  select(., SampleID, fitness, preadapt_pop, evolution) %>%
  rename(., Var1 = SampleID)
d_meta_cont <- sample_data(ps_prop) %>% 
  data.frame() %>%
  filter(., treatment != 'individual_clone') %>%
  select(., SampleID, treatment) %>%
  rename(., Var2 = SampleID)

d_wunifrac <- filter(d_wunifrac, Var1 != Var2) %>%
  filter(., Var1 %in% d_meta_ind_clone$Var1) %>%
  filter(., Var2 %in% d_meta_cont$Var2)
d_wunifrac <- merge(d_wunifrac, d_meta_cont, by = 'Var2') %>% merge(., d_meta_ind_clone, by = 'Var1') %>%
  rename(., wunifrac = Freq)

ggplot(d_wunifrac, aes(fitness, wunifrac, col = evolution)) +
  geom_point() +
  facet_wrap(~ treatment)
# NOPE - definitely not

d_fitness <- filter(d_wunifrac, treatment == 'lacz_ancestor') %>%
  distinct(., Var1, .keep_all = TRUE)

ggplot(d_fitness, aes(evolution, fitness)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.1, height = 0.1), size = 3) +
  theme_bw(base_size = 14) +
  ggtitle('Effect of preadaptation history on individual clone fitness')

ggsave(file.path(path_fig, 'ind_clone_fitness.png'), last_plot(), height = 5, width = 6)

mod <- lm(fitness ~ evolution, d_fitness)
mod1 <- lmer(fitness ~ evolution + (1|preadapt_pop), d_fitness)
mod2 <- lmer(fitness ~ 1 + (1|preadapt_pop), d_fitness)

mod1 <- lmer(Freq ~ fitness*treatment + (1|Var1), d_wunifrac, na.action = na.fail, REML = FALSE)
MuMIn::dredge(mod1)
mod2 <- lmer(Freq ~ fitness + treatment + (1|Var1), d_wunifrac, na.action = na.fail, REML = FALSE)
mod3 <- lmer(Freq ~ treatment + (1|Var1), d_wunifrac, na.action = na.fail, REML = FALSE)


# analysis of just individual clones, wild type (lacz) vs pre-adapted ####

# filter samples that are not individual clones ####
# specifically negative control, nmc_T0 & wt_ancestor (dont care as only lacz were in the communities)
# specifically wt ancestor and nmc_t0
to_keep <- filter(meta_new, treatment %in% c('individual_clone', 'lacz_ancestor'))
ps2 <- prune_samples(to_keep$SampleID, ps)

# transform counts to relative abundances for ordination ####
ps_prop <- transform_sample_counts(ps2, function(x){x / sum(x)})

# Weighted Unifrac ####
# weighted Unifrac distance
# to change the ordination - see https://joey711.github.io/phyloseq/plot_ordination-examples.html and https://joey711.github.io/phyloseq/distance.html

# do weighted unifrac
ord_wUni <- ordinate(ps_prop, method = 'MDS', distance = 'wunifrac')

evals <- ord_wUni$values$Eigenvalues

# plot
plot_ordination(ps_prop, ord_wUni, color = "treatment") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~ treatment)
# plot all on one panel by removing the facet_wrap command

# get the distance matrix out of the data
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp,
                 treatment = as.factor(treatment)) %>%
  column_to_rownames(., 'SampleID')

# run an Adonis test
mod_treatment <-  vegan::adonis(ps_wunifrac ~ treatment, data = d_samp, n_perm = 9999)

# overwrite metadata to allow plotting of nclones_fac
sample_data(ps_prop) <- sample_data(d_samp)

# plot ordination
plot_ordination(ps_prop, ord_wUni, color = 'treatment') +
  geom_point(size = 2) +
  stat_ellipse(aes(fill = treatment, group = treatment), geom = 'polygon', type = "t", alpha = 0.05) +
  theme_bw(base_size = 10, base_family = 'Helvetica') +
  ylab('PCoA2 [19.6%]') +
  xlab('PCoA1 [44.7%]')

# beta-diversity analysis - look at homogeneity of variances
mod1_dispers <- betadisper(ps_wunifrac, d_samp$treatment)

# plot of model
plot(mod1_dispers)
boxplot(mod1_dispers)

# anova
anova(mod1_dispers)

# Permutation test for F
pmod <- permutest(mod1_dispers, pairwise = TRUE)

# Tukey's Honest Significant Differences
T_HSD <- TukeyHSD(mod1_dispers)

# plot distances to centroid - try and play with alpha

# get betadisper data ####
betadisper_dat <- get_betadisper_data(mod1_dispers)

# do some transformations on the data
betadisper_dat$eigenvalue <- mutate(betadisper_dat$eigenvalue, percent = eig/sum(eig)*100)

# add convex hull points ####
# this could be put in a function
betadisper_dat$chull <- group_by(betadisper_dat$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(betadisper_dat$centroids, group, PCoA1, PCoA2), select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
betadisper_dat$eigenvector$distances <- betadisper_dat$distances$distances

# plot PCoA
p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 5) +
  geom_point(aes(PCoA1, PCoA2, col = group, alpha = 1 - distances), betadisper_dat$eigenvector, size = 0.75) +
  #geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull, alpha = 0.1) +
  #geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = group, alpha = 1 - distances), betadisper_lines) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  stat_ellipse(aes(PCoA1, PCoA2, col = group, group = group), type = "t", betadisper_dat$eigenvector) +
  ylab('PCoA Axis 2 [18.37%]') +
  xlab('PCoA Axis 1 [44.31%]') +
  theme(legend.position = 'bottom') +
  coord_fixed(sqrt(betadisper_dat$eigenvalue$percent[2]/betadisper_dat$eigenvalue$percent[1])) +
  ggtitle('PCoA plot with closer points emphasised') +
  scale_alpha(range = c(0.0001, 0.75), guide = FALSE)

# plot distances from centroid
p2 <- ggplot(betadisper_dat$distances, aes(group, distances, fill = group, col = group)) +
  geom_boxplot(aes(fill = group, col = group), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(group, distances, col = group), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.3)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('') +
  ggtitle('Distance from centroid')

p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(2.5, 1))

ggsave(file.path(path_fig, 'wt_vs_preadapt.pdf'), p3, height = 5, width = 10)

# filter samples that are not individual clones ####
# specifically negative control, nmc_T0 & wt_ancestor (dont care as only lacz were in the communities)
# specifically wt ancestor and nmc_t0
to_keep <- filter(meta_new, treatment %in% c('evolved_with_community', 'evolved_without_community', 'lacz_ancestor'))
ps2 <- prune_samples(to_keep$SampleID, ps)

# transform counts to relative abundances for ordination ####
ps_prop <- transform_sample_counts(ps2, function(x){x / sum(x)})

# Weighted Unifrac ####
# weighted Unifrac distance
# to change the ordination - see https://joey711.github.io/phyloseq/plot_ordination-examples.html and https://joey711.github.io/phyloseq/distance.html

# do weighted unifrac
ord_wUni <- ordinate(ps_prop, method = 'MDS', distance = 'wunifrac')

evals <- ord_wUni$values$Eigenvalues

# plot
plot_ordination(ps_prop, ord_wUni, color = "n_clones") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~ treatment)
# plot all on one panel by removing the facet_wrap command

# get the distance matrix out of the data
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp,
                 n_clones_fac = as.factor(n_clones),
                 treatment = as.factor(treatment)) %>%
  column_to_rownames(., 'SampleID')

# run an Adonis test
mod_nclones <-  vegan::adonis(ps_wunifrac ~ n_clones_fac, data = d_samp, n_perm = 9999)
mod_treat <-  vegan::adonis(ps_wunifrac ~ treatment, data = d_samp, n_perm = 9999)

# run a multiple comparison to see which treatments are different
mult_comp <- mctoolsr::calc_pairwise_permanovas(ps_wunifrac, d_samp, 'treatment', n_perm = 9999)
# loads of comparisons. Significant differences will be determined by p value correction.

# overwrite metadata to allow plotting of nclones_fac
sample_data(ps_prop) <- sample_data(d_samp)

# beta-diversity analysis - look at homogeneity of variances
mod1_dispers <- betadisper(ps_wunifrac, d_samp$n_clones_fac)

# plot of model
plot(mod1_dispers)
boxplot(mod1_dispers)

# anova
anova(mod1_dispers)

# Permutation test for F
pmod <- permutest(mod1_dispers, pairwise = TRUE)

# Tukey's Honest Significant Differences
T_HSD <- TukeyHSD(mod1_dispers)

# plot distances to centroid - try and play with alpha

# get betadisper data ####
betadisper_dat <- get_betadisper_data(mod1_dispers)

# do some transformations on the data
betadisper_dat$eigenvalue <- mutate(betadisper_dat$eigenvalue, percent = eig/sum(eig)*100)

# add convex hull points ####
# this could be put in a function
betadisper_dat$chull <- group_by(betadisper_dat$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(betadisper_dat$centroids, group, PCoA1, PCoA2), select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
betadisper_dat$eigenvector$distances <- betadisper_dat$distances$distances

# plot PCoA
p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 5) +
  geom_point(aes(PCoA1, PCoA2, col = group, alpha = 1 - distances), betadisper_dat$eigenvector, size = 0.75) +
  #geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull, alpha = 0.1) +
  #geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = group, alpha = 1 - distances), betadisper_lines) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  stat_ellipse(aes(PCoA1, PCoA2, col = group, group = group), type = "t", betadisper_dat$eigenvector) +
  ylab('PCoA Axis 2 [19.03%]') +
  xlab('PCoA Axis 1 [46.78%]') +
  theme(legend.position = 'bottom') +
  coord_fixed(sqrt(betadisper_dat$eigenvalue$percent[2]/betadisper_dat$eigenvalue$percent[1])) +
  ggtitle('PCoA of lacz ancestor vs 24 clone preadapted mixtures') +
  scale_alpha(range = c(0.0001, 0.75), guide = FALSE)

# plot distances from centroid
p2 <- ggplot(betadisper_dat$distances, aes(group, distances, fill = group, col = group)) +
  geom_boxplot(aes(fill = group, col = group), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(group, distances, col = group), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.3)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('') +
  ggtitle('Distance from centroid')

p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(2.5, 1))

ggsave(file.path(path_fig, 'wt_vs_populations.pdf'), p3, height = 5, width = 10)



