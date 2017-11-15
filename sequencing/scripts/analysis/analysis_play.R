# play with analyses of sequencing data ####
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
library(DESeq2)

# set seed
set.seed(42)

# figure path
path_fig <- 'sequencing/plots'

# load data - latest run which we are happy with ####
ps <- readRDS('sequencing/data/output/20171024_17:18/20171024_17:18_ps.rds')

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
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~ treatment)

plot_ordination(ps_prop, ord_wUni, type = 'taxa', color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~ Phylum)

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

################################
# create taxonomy summaries ####

# summarise phyloseq object at the Phylum level
phylum_glom <- tax_glom(ps, taxrank = "Phylum" )

# filter for the 10 most common phyla

# function for ease of use
# code stolen from phyloseq website
get_top_taxa <- function(ps, tax_rank, to_keep){
  temp <- tapply(phyloseq::taxa_sums(ps), phyloseq::tax_table(ps)[, tax_rank], sum, na.rm = TRUE)
  temp2 <-  names(sort(temp, TRUE))[1:to_keep]
  return(temp2)
}

# filter
phylum_glom_filt <- prune_taxa((tax_table(phylum_glom)[, "Phylum"] %in% get_top_taxa(phylum_glom, 'Phylum', 10)), phylum_glom)

# convert counts to proportions
phylum_prop <- transform_sample_counts(phylum_glom, function(x){x / sum(x)})
phylum_prop_filt <- transform_sample_counts(phylum_glom_filt, function(x){x / sum(x)})

# plot bar plot
plot_bar(phylum_prop_filt, fill = "Phylum") +
  facet_wrap(~ treatment, scale = 'free_x')

# try to create a prettier bar plot

# get data
d_glom <- psmelt(phylum_glom_filt)

# group by treatments
d_glom_group <- group_by(d_glom, treatment, preadapt_pop, evolution) %>%
  do(., data.frame(prop = .$Abundance/sum(.$Abundance), Phylum = .$Phylum)) %>%
  ungroup()

ggplot(d_glom_group, aes(treatment, prop, fill = Phylum)) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ preadapt_pop, scale = 'free_x')

# re ordinate the data!!! ####
#####################

# double principal component analysis
ord_dpcoa <- ordinate(ps_prop, method = 'DPCoA')

evals <- ord_dpcoa$eig

p_taxa_dpcoa <- plot_ordination(ps_prop, ord_dpcoa, type = 'taxa', color = "Phylum") +
  #coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('DPCoA plot of phylum position') 

p_samp_dpcoa <- plot_ordination(phylum_prop_filt, ord_dpcoa, type = 'samples', color = "treatment") +
  #coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('DPCoA plot of ')

grid.arrange(p_taxa_dpcoa, p_samp_dpcoa)

# Weighted Unifrac ####
ord_wUni <- ordinate(phylum_prop, method = 'MDS', distance = 'wunifrac')

evals <- ord_wUni$values$Eigenvalues

plot_ordination(phylum_prop, ord_wUni, color = "treatment") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~treatment)

plot_ordination(four_clones, ord_wUni, color = "preadapt_pop", shape = 'evolution') +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap(~ evolution)




# individual clones only ####
ind_clonesID <- filter(meta_new, treatment == 'individual_clone')
ind_clones <- prune_samples(ind_clonesID$SampleID, phylum_prop)

ord_wUni <- ordinate(ind_clones, method = 'MDS', distance = 'wunifrac')

evals <- ord_wUni$values$Eigenvalues

plot_ordination(ind_clones, ord_wUni, color = "preadapt_pop", shape = 'evolution') +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap(~ evolution)

# 4 clones only 
# individual clones only ####
four_clonesID <- filter(meta_new, treatment %in% c('4_related_clones', '4_unrelated_clones')) %>%
  mutate(., evolution = ifelse(is.na(evolution), 'both', evolution))
four_clones <- prune_samples(four_clonesID$SampleID, phylum_prop)
sample_data(four_clones)$evolution <- ifelse(is.na(sample_data(four_clones)$evolution), 'both', sample_data(four_clones)$evolution)
sample_data(four_clones)$preadapt_pop <- ifelse(is.na(sample_data(four_clones)$preadapt_pop), 'many', sample_data(four_clones)$preadapt_pop)

ord_wUni <- ordinate(four_clones, method = 'MDS', distance = 'wunifrac')

evals <- ord_wUni$values$Eigenvalues

plot_ordination(four_clones, ord_wUni, color = "preadapt_pop", shape = 'evolution') +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Unifrac distances') +
  facet_wrap(~ evolution)



plot_net(ps, maxdist = 0.4, point_label = "treatment", col = 'evolution')

# quick play with Deseq ####
d_deseq <- phyloseq_to_deseq2(ps, ~ treatment)
d_deseq <- DESeq(d_deseq, test = 'Wald', fitType = 'parametric')

# results table
res = results(d_deseq, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

# plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))