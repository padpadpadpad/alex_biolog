# analysis trying to predict  ####
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
library(ggvegan)
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
ps <- readRDS('sequencing/data/output/20171024_17:18/20171024_17:18_ps.rds')

# replace metadata with new metadata
# when wanting to add columns to metadata, it is better to edit metadata_creation and overwrite the metadata file as then it can be overwritten in all future files
meta_new <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE)
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 30,000. Woof.

# not going to rarefy those samples yet

# DIVERSITY ANALYSIS ####
# Look at the effect of number of clones on community composition (as clustering)

# filter some samples with NA for n_clones ####
# specifically negative control, nmc_T0 & wt_ancestor (dont care as only lacz were in the communities)
# specifically wt ancestor and nmc_t0
to_keep <- filter(meta_new, treatment %in% c('individual_clone', '4_related_clones', '4_unrelated_clones'))
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
plot_ordination(ps_prop, ord_wUni, col = 'n_clones') +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point() +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~ preadapt_pop)
# plot all on one panel by removing the facet_wrap command



# save plot, other ways are available
ggsave(file.path(path_fig, 'ordination.pdf'), last_plot(), height = 6, width = 12)

# get the distance matrix out of the data
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# Need to do a PCoA to get the weighted unifrac distances into Euclidean space. capscale() should help with this
pcoa <- capscale(ps_wunifrac ~ 1)


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
  ylab('PCoA2 [18.4%]') +
  xlab('PCoA1 [44.2%]')
