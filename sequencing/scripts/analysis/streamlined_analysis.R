# streamlined script for analysis thus far ####
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

# set seed
set.seed(42)

# figure path
path_fig <- 'sequencing/plots'

# load data - latest run which we are happy with ####
# these files need to be there
ps <- readRDS('sequencing/data/output/20171024_17:18/20171024_17:18_ps.rds')

# replace metadata with new metadata
meta_new <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE)
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 30,000. Woof.

# not going to rarefy those samples yet

# can plot rarefaction curves
# code from https://www.fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
# written by Gavin Simpson himself (an author of vegan)

# check rarefaction curves ####
ps_otu_table <- data.frame(otu_table(ps))
raremax <- min(rowSums(ps_otu_table))
col <- c('black', 'blue', 'yellow', 'red', 'orange', 'grey', 'hotpink', 'purple', 'green')
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
out <- with(pars,
            rarecurve(ps_otu_table, step = 1000, sample = raremax, col = col,
                      lty = lty, label = FALSE))

# having done similar analyses on the rarefied samples they give very similar results. Therefore will use the unrarefied samples for the ordinations

# filter some samples ####
# specifically wt ancestor and nmc_t0
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0', 'wt_ancestor'))
ps <- prune_samples(to_keep$SampleID, ps)

# transform counts to relative abundances for ordination ####
ps_prop <- transform_sample_counts(ps, function(x){x / sum(x)})

# weighted Unifrac distance
# to change the ordination - see https://joey711.github.io/phyloseq/plot_ordination-examples.html and https://joey711.github.io/phyloseq/distance.html

# 1. Weighted Unifrac ####

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
d_samp$treatment <- as.factor(d_samp$treatment)

# run an Adonis test
mod1 <- vegan::adonis(ps_wunifrac ~ treatment, data = d_samp)

# run a multiple comparison to see which treatments are different
mult_comp <- mctoolsr::calc_pairwise_permanovas(ps_wunifrac, d_samp, 'treatment', n_perm = 9999)
# loads of comparisons. Significant differences will be determined by p value correction.

################################
# create taxonomy summaries ####

# summarise phyloseq object at the Phylum level
tax_group <- tax_glom(ps, taxrank = "Genus")

# filter for the 10 most common genus ####
tax_group_filt <- prune_taxa((tax_table(tax_group)[, "Genus"] %in% get_top_taxa(tax_group, 'Genus', 10)), tax_group)

# convert counts to proportions
tax_prop <- transform_sample_counts(tax_group, function(x){x / sum(x)})
tax_prop_filt <- transform_sample_counts(tax_group_filt, function(x){x / sum(x)})

# plot bar plot
plot_bar(tax_prop_filt, fill = "Genus") +
  facet_wrap(~ treatment, scale = 'free_x')

# save plot, other ways are available
ggsave(file.path(path_fig, 'plot_bar1.pdf'), last_plot(), height = 7, width = 12)

# try to create a prettier bar plot

# get data
d_glom <- psmelt(tax_prop_filt)

# group by treatments
d_glom_group <- group_by(d_glom, treatment, preadapt_pop, evolution) %>%
  do(., data.frame(prop = .$Abundance/sum(.$Abundance), Genus = .$Genus)) %>%
  ungroup()

# plot these
ggplot(d_glom_group, aes(treatment, prop, fill = Genus, col = Genus)) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('Relative abundance of 10 most abundant genus groups')

ggsave(file.path(path_fig, 'plot_bar2.pdf'), last_plot(), height = 5, width = 8)

# re ordinate the data to see if aggregating each phylum made a difference ####

# Weighted Unifrac ####
ord_wUni2 <- ordinate(tax_prop, method = 'MDS', distance = 'wunifrac')

evals2 <- ord_wUni2$values$Eigenvalues

plot_ordination(tax_prop, ord_wUni2, color = "treatment") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~treatment)

# plot original again
plot_ordination(ps_prop, ord_wUni, color = "treatment") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~ treatment)

# does change things - do not think this is a huge problem right now

# do beta diversity analysis ####
# Looking at variance across groups
# assume all samples are independent

# use raw proportion data

# can only do one factor would use every combination of id
# have negative eigenvalues, need to correct for these...?
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

# no differences...

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

# Now the dataframes are all ready to be completely customisable in ggplot
# plot betadispersion plot
ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 4) +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$eigenvector) +
  geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull ) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = group), betadisper_lines) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('PCoA Axis 2 [17.6%]') +
  xlab('PCoA Axis 2 [45.6%]') +
  theme(legend.position = c(0.9, 0.1)) +
  coord_fixed(sqrt(betadisper_dat$eigenvalue$percent[2]/betadisper_dat$eigenvalue$percent[1])) +
  facet_wrap(~ group) +
  ggtitle('PCoA plot with centroids per treatment')

ggsave(file.path(path_fig, 'PCoA_plot.pdf'), last_plot(), height = 6, width = 12)

# plot distances from centroid
ggplot(betadisper_dat$distances, aes(group, distances, fill = group, col = group)) +
  geom_boxplot(aes(fill = group, col = group), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(group, distances, col = group), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.2)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('Distance to centroid') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('') +
  ggtitle('Distance from centroid')

ggsave(file.path(path_fig, 'dist_centroid.pdf'), last_plot(), height = 5, width = 8)

# plot as joy plot
ggplot(betadisper_dat$distances, aes(distances, group, fill = group, col = group)) +
  geom_density_ridges() +
  xlab('distance from centroid') +
  ggtitle('Distance from centroid')

ggsave(file.path(path_fig, 'dist_centroid2.pdf'), last_plot(), height = 6, width = 10)

