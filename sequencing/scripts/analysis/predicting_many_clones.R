# analysis trying to predict  ####
rm(list = ls())

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vegan)
library(tibble)
library(patchwork) # devtools::install_github('thomasp85/patchwork')
library(ggvegan)
# if not installed, install mctoolsr run devtools::install_github('leffj/mctoolsr')

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
ps_bray <- phyloseq::distance(ps2, method = 'bray')

# Need to do a PCoA to get the weighted unifrac distances into Euclidean space. capscale() should help with this
pcoa <- capscale(ps_bray ~ 1)
plot(pcoa)

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp_ind <- filter(d_samp, treatment == 'individual_clone')
d_samp_4 <- filter(d_samp, treatment != 'individual_clone')

# create column for clone
d_samp_ind <- mutate(d_samp_ind, clone = paste('C', sample_name, sep = '_')) %>%
  mutate(., rep = group_indices(., preadapt_pop),
         evol = 'related')

# try and get the data of the pcoa out
d_pcoa <- fortify(pcoa) %>%
  janitor::clean_names() %>%
  rename(., SampleID = label) %>%
  merge(., select(d_samp, SampleID, sample_name))

# metaMDS - nmds
nmds <- metaMDS(ps_bray)
stressplot(nmds)

d_nmds <- data.frame(scores(nmds)) %>%
  janitor::clean_names() %>%
  tibble::rownames_to_column(., 'SampleID') %>%
  merge(., select(d_samp, SampleID, sample_name))

# choose ordination method to use ###
d_ord <- d_nmds

# filter for just the clones
d_ord_ind <- filter(d_ord, sample_name %in% as.character(1:48)) %>%
  mutate(., clone = paste('C', sample_name, sep = '_'))

# get every combination used in the experiment

# read in mixed population data
d_samp_4_ind <- read.csv('sequencing/data/Mixed_population_data.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  mutate(., clone = paste('C', clone, sep = '_'),
         rep = group_indices(., preadaptation_population),
         preadapt_pop = preadaptation_population,
         evol = 'unrelated') %>%
  select(., -preadaptation_population)

# bind with d_samp_ind
d_ind_combs <- rbind(select(d_samp_ind, sample_name, clone, rep, preadapt_pop, evol),
                     d_samp_4_ind)
# This is good

# merge PCoA axes with clone
d_ord_ind <- merge(d_ind_combs, select(d_ord_ind, c(clone, contains('mds'))), by = 'clone', all.x = TRUE)
# This seems to have worked

# add a unique ID for each population
d_ord_ind <- mutate(d_ord_ind, rep2 = case_when(sample_name %in% as.character(1:48) ~ preadapt_pop,
                                                  TRUE ~ sample_name))

# calculate centroid of each population - using the "mediod" as this is the method used in beta_disper, could just as easily get the mean

cols <- Hmisc::Cs(mds1, mds2, mds3, mds4, mds5, mds6, nmds1, nmds2)

d_ord_centroids <- group_by(d_ord_ind, rep2) %>%
  summarise_at(., vars(contains('mds')), median) %>%
  ungroup()

# filter just the 4 clones of the pcoa
d_ord_4 <- filter(d_ord, ! SampleID %in% as.character(1:48)) %>%
  merge(., select(d_samp_4, SampleID, preadapt_pop), by = 'SampleID') %>%
  mutate(., rep2 = case_when(is.na(preadapt_pop) ~ sample_name,
                             TRUE ~ preadapt_pop)) %>%
  select(., -preadapt_pop)

# first plot
ggplot(d_ord_ind, aes(nmds1, nmds2)) +
  geom_point() +
  geom_point(data = d_ord_centroids, col = 'red', size = 2.5) +
  geom_point(data = d_ord_4, col = 'blue', shape = 2) +
  facet_wrap(~ rep2)

# calculate average distance between centroid and corresponding sample ####
 
# for WT-B6 there is no sample of the four clones together
# delete WT-B6 from list

id <- unique(d_ord_ind$rep2)
id <- id[id != 'WT-B6']

# distance is calculated as
# sqrt(sum((x_i - y_i)^2)) - Euclidean distance

d_ord_centroids <- mutate(d_ord_centroids, samp = 'centroid')
d_ord_4 <- mutate(d_ord_4, samp = '4_clones')

common_cols <- intersect(colnames(d_ord_centroids), colnames(d_ord_4))

# bind centroid and common columns together
d_dist <- bind_rows(select(d_ord_4, common_cols), select(d_ord_centroids, common_cols))

# get rid of data which cannot be used for comparison
d_dist <- filter(d_dist, rep2 != 'WT-B6')

d_dist_sum <- gather(d_dist, 'eigenvector', 'value', contains('mds')) %>%
  spread(., samp, value) %>%
  group_by(., rep2, eigenvector) %>%
  summarise(., dist = (`4_clones` - centroid)^2) %>%
  ungroup() %>%
  group_by(., rep2) %>%
  summarise(., distance = sqrt(sum(dist))) %>%
  ungroup()

# average distance of actual samples and replicates
actual_dist = mean(d_dist_sum$distance, na.rm = TRUE)  

# permute this distance many many times ####

# number of permutations
n_perm = 9999

# create empty vector
perm_dist <- rep(NA, times = n_perm)

# set progress bar
pb <- progress::progress_bar$new(total = n_perm, clear = FALSE)

# run for loop for shuffling
for(i in 1:n_perm){
  
  # progress bar
  pb$tick()
  # shuffle up the samples, randomly within centroids and 4 clones respectively (I think)
  # while keeping all the mds' the same
  d_dist_temp <- group_by(d_dist, samp) %>%
    mutate(., rep2 = sample(rep2, replace = FALSE)) %>%
    ungroup()
  
  # run the code to calculate Euclidean distance
  d_dist_sum_temp <- gather(d_dist_temp, 'eigenvector', 'value', contains('mds')) %>%
    spread(., samp, value) %>%
    group_by(., rep2, eigenvector) %>%
    summarise(., dist = (`4_clones` - centroid)^2) %>%
    ungroup() %>%
    group_by(., rep2) %>%
    summarise(., distance = sqrt(sum(dist))) %>%
    ungroup()
  
  # save out the mean distance between centroid and 4 clone performance
  perm_dist[i] <- mean(d_dist_sum_temp$distance, na.rm = TRUE)
  
  # remove unwanted objects
  rm(list = c('d_dist_temp', 'd_dist_sum_temp'))
}

# plot this ####

# second plot
ggplot(d_ord_ind, aes(nmds1, nmds2)) +
  geom_point(size = 0.5) +
  geom_point(aes(shape = samp), data = d_dist, col = 'red', size = 2.5) +
  geom_line(aes(group = rep2), data = d_dist, col = 'red', alpha = 0.5) +
  theme(legend.position = c(0.9, 0.9)) +
  coord_equal()

# plot with facets
p1 <- ggplot(d_ord_ind, aes(nmds1, nmds2)) +
  geom_point(size = 0.5) +
  geom_point(aes(shape = samp), data = d_dist, col = 'red', size = 2.5) +
  geom_line(aes(group = rep2), data = d_dist, col = 'red', alpha = 0.5) +
  theme(legend.position = 'none') +
  coord_equal()+
  facet_wrap(~rep2) +
  ggtitle('Position of each replicate in euclidean space based on an NMDS',
          subtitle = 'Red circles are 4 clone responses, triangles are centroid positions') +
  theme_bw(base_size = 8)

d_perm_dist = data.frame(n_perm = 1:n_perm, dist = perm_dist)

p2 <- ggplot(d_perm_dist, aes(dist)) +
  geom_histogram(fill = 'white', col = 'black') +
  geom_vline(aes(xintercept = actual_dist), col = 'red') +
  theme_grey(base_size = 8) +
  xlab('Euclidean distance between two points') +
  ggtitle('Permuted euclidean distance between centroid of 4 individual clones and actual effect of 4 clones',
          subtitle = 'Red line is average distance of actual measured data')

# calculate permuted p-value ####
# How many distances are less than the true distance
# p = sum(perm_dist - actual_dist < 0) /n_perm + 1
sum(perm_dist - actual_dist < 0) /(n_perm + 1)

# with NMDS, the pvalue = 0.1149
p3 <- p1 + p2
p3

ggsave(file.path(path_fig, 'predict_multiclones.png'), p3, height = 5, width = 14)
