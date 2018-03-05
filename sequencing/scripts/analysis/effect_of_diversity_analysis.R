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
# if not installed, install mctoolsr run devtools::install_github('leffj/mctoolsr')

# functions ####
# code stolen from phyloseq website
get_top_taxa <- function(ps, tax_rank, to_keep){
  temp <- tapply(phyloseq::taxa_sums(ps), phyloseq::tax_table(ps)[, tax_rank], sum, na.rm = TRUE)
  temp2 <-  names(sort(temp, TRUE))[1:to_keep]
  return(temp2)
}

# get betadisper dataframes
# have written functions to grab the necessary data from the betadisper object

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

# set seed
set.seed(42)

# figure path
path_fig <- 'sequencing/plots'

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

# not going to rarefy those samples yet

# DIVERSITY ANALYSIS ####
# Look at the effect of number of clones on community composition (as clustering)

# filter some samples with NA for n_clones ####
# specifically negative control, nmc_T0 & wt_ancestor (dont care as only lacz were in the communities)
# specifically wt ancestor and nmc_t0
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0', 'negative_control', 'wt_ancestor'))
ps2 <- prune_samples(to_keep$SampleID, ps)

# remove Pseudomonas reads from the analysis ####
# have a look at number of genus 
table(tax_table(ps2)[, "Genus"], exclude = NULL)
# have a look at Pseudomonas more closely
d_pseu <- transform_sample_counts(ps2, function(x){x / sum(x)}) %>%
  subset_taxa(., Genus == 'Pseudomonas') %>%
  psmelt() %>%
  group_by(SampleID) %>%
  summarise(prop = sum(Abundance))
# no species so delete all pseudomonas samples from data

# remove pseudomonas
ps2 <- subset_taxa(ps2, Genus != 'Pseudomonas')

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

# save plot, other ways are available
ggsave(file.path(path_fig, 'ordination.pdf'), last_plot(), height = 6, width = 12)

# get the distance matrix out of the data
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac)) %>%
  column_to_rownames(., 'SampleID')

# run an Adonis test
mod_nclone <- vegan::adonis(ps_wunifrac ~ n_clones, data = d_samp, n_perm = 9999)
mod_nclonefac <-  vegan::adonis(ps_wunifrac ~ nclones_fac, data = d_samp, n_perm = 9999)

# run a multiple comparison to see which treatments are different
mult_comp <- mctoolsr::calc_pairwise_permanovas(ps_wunifrac, d_samp, 'nclones_fac', n_perm = 9999)
# loads of comparisons. Significant differences will be determined by p value correction.

# overwrite metadata to allow plotting of nclones_fac
sample_data(ps_prop) <- sample_data(d_samp)

# plot ordination
plot_ordination(ps_prop, ord_wUni, color = "nclones_fac", shape = 'treatment') +
  geom_point(size = 2) +
  scale_shape_discrete('Treatment') +
  stat_ellipse(aes(fill = nclones_fac, group = nclones_fac), geom = 'polygon', type = "t", alpha = 0.05) +
  theme_bw(base_size = 10, base_family = 'Helvetica') +
  ylab('PCoA2 [23%]') +
  xlab('PCoA1 [45.9%]')

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
p2 <- ggplot(betadisper_dat$distances, aes(forcats::fct_relevel(group, 'C_24', after = 2), distances, fill = group, col = group)) +
  geom_boxplot(aes(fill = group, col = group), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(group, distances, col = group), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.3)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('') +
  ggtitle('Distance from centroid')

p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(2.5, 1))

ggsave(file.path(path_fig, 'PCoA_plot_diversity.pdf'), p3, height = 5, width = 10)

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



