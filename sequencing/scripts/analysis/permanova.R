# PERMANOVA and ordination ####
rm(list = ls())

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vegan)
library(gridExtra)

# figure path
path_fig <- 'sequencing/plots'

# load data - latest run which we are happy with ####
ps <- readRDS('sequencing/data/output/20171024_17:18/20171024_17:18_ps.rds')

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

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps, function(x){x / sum(x)})

# 1. unweighted Unifrac analysis ####
ord_unUni <- ordinate(ps_prop, method = 'MDS')

evals <- ord_unUni$values$Eigenvalues

p_unUni <- plot_ordination(ps_prop, ord_unUni, color = "treatment") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on Unifrac distances')

# calculate distance matrix
ps_unifrac <- phyloseq::distance(ps_prop, method = 'unifrac')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps))

# run an Adonis test
mod1 <- vegan::adonis(ps_unifrac ~ treatment, data = d_samp)

#####################
# 2. Weighted Unifrac ####
ord_wUni <- ordinate(ps_prop, method = 'MDS', distance = 'wunifrac')

evals <- ord_wUni$values$Eigenvalues

p_wUni <- plot_ordination(ps_prop, ord_wUni, color = "treatment") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances')

ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))

# run an Adonis test
vegan::adonis(ps_wunifrac ~ treatment, data = d_samp)

#####################
# 3. Bray Curtis ####
ord_Bray <- ordinate(ps_prop, method = 'NMDS', distance = 'bray', trymax = 100, k = 2)

evals <- ord_Bray$values$Eigenvalues

p_Bray <- plot_ordination(ps_prop, ord_Bray, color = "treatment") +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('NMDS of bray curtis similarity')

ps_bray <- phyloseq::distance(ps_prop, method = 'bray')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))

# run an Adonis test
vegan::adonis(ps_bray ~ treatment, data = d_samp)

# plot all three plots
plot_all <- grid.arrange(p_unUni + theme(legend.position = 'none'), 
             p_wUni + theme(legend.position = 'none'), 
             p_Bray, 
             ncol = 3)

ggsave(file.path(path_fig, 'ordination.pdf'), plot_all, height = 7, width = 20)
