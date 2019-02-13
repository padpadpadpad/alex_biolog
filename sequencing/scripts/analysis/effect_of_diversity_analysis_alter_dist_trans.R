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
library(patchwork)
library(lme4)
library(mctoolsr)

#--------------#
# functions ####
#--------------#

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

# convert distance matrix to dataframe
dist_2_df <- function(dist_ob){
  m <- as.matrix(dist_ob) # coerce dist object to a matrix
  xy <- t(combn(colnames(m), 2))
  return(data.frame(xy, dist=m[xy], stringsAsFactors = FALSE))
}

#-------------------------------------#
# setup workspace and load in data ####
#-------------------------------------#

# set seed
set.seed(42)

# figure path
path_fig <- 'plots/validation'

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

# wrap this code in a function

f_analysis <- function(phyloseq_object, metric = 'bray', abund = 'raw'){
  
  # create file
  file.create(file.path('results', paste(abund, metric, 'diversity_output.txt', sep = '_')))
  file <- file(file.path('results', paste(abund, metric, 'diversity_output.txt', sep = '_')))
  
  # sink all the output to this file (will put all console output into this file)
  sink(file, type = 'output')
  
  # print out the start line
  cat(paste('###########################################\n# This is the output file from the sequencing analysis for the Pseudomonas paper when the transformation used was:', abund, 'and the distance metric used was:', metric, '\n###########################################\n', sep = ' '), file = file, append = TRUE)
  
  # do transformation
  ps_transform <- phyloseq_object
  
  # supply methods for transformations
  if(abund == 'prop'){ps_transform <- transform_sample_counts(ps_transform, function(x){x / sum(x)})}
  if(abund == 'hellinger'){ps_transform <- transform_sample_counts(ps_transform, function(x){sqrt(x / sum(x))})}
  if(abund == 'log'){ps_transform <- transform_sample_counts(ps_transform, function(x){log(x + 1)})}
  if(abund == 'sqrt'){ps_transform <- transform_sample_counts(ps_transform, function(x){sqrt(x)})}
  
  
  # -------------------------------------------------------------------------------#
  # Analysis of the preadaptation with and without nmc on community composition ####
  #--------------------------------------------------------------------------------#
  
  # perfectly replicated for 1, 4 and 24 clones
  # do these separately and then pvalue correct
  
  # 1. individual clones ####
  
  # samples to keep - just individual clones
  to_keep <- filter(meta_new, treatment %in% c('individual_clone'))
  ps_sub <- prune_samples(to_keep$SampleID, ps_transform)
  
  # make a data frame of the sample data
  d_samp <- data.frame(sample_data(ps_sub))
  d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                   nclones_fac = as.factor(nclones_fac),
                   evol_fac = as.factor(evolution),
                   preadapt_pop = as.factor(preadapt_pop),
                   pop2 = preadapt_pop) %>%
    column_to_rownames(., 'SampleID')
  levels(d_samp$pop2) <- letters[1:length(levels(d_samp$pop2))]
  
  # calculate distance matrix
  ps_dist <- phyloseq::distance(ps_sub, method = metric)
  
  # calculate average distance of each clone to its sympatyric populations and allopatric populations across reps
  d_dist <- dist_2_df(ps_dist) %>%
    merge(., select(d_samp, X1_evol = evolution, X1_preadapt_pop = preadapt_pop) %>% rownames_to_column(., var = 'X1'), by = 'X1') %>%
    merge(., select(d_samp, X2_evol = evolution, X2_preadapt_pop = preadapt_pop) %>% rownames_to_column(., var = 'X2'), by = 'X2') %>%
    filter(., X1_evol == X2_evol) %>%
    mutate(same_rep = ifelse(X1_preadapt_pop == X2_preadapt_pop, 'Y', 'N')) %>%
    gather(., 'random', 'clone', c(X1, X2)) %>%
    group_by(clone, same_rep, X1_evol) %>%
    summarise(., mean = mean(dist)) %>%
    ungroup()
  
  
  labels <- c(with_community = '(a) pre-adapted with nmc', without_community = '(b) pre-adapted without nmc')
  
  p1 <- ggplot(d_dist, aes(same_rep, mean)) +
    MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
    geom_point(fill = 'white', shape = 21, size = 3, position = position_jitter(width = 0.1)) +
    facet_wrap(~ X1_evol, labeller = labeller(X1_evol = labels)) +
    theme_bw(base_size = 16) +
    ylab(paste(metric, 'distance', sep = ' ')) +
    xlab('') +
    scale_x_discrete(labels = c('allopatric pair', 'sympatric pair')) +
    theme(strip.background = element_blank(),
          strip.text = element_text(hjust = 0)) +
    ggtitle('Average distance between allopatric and sympatric clones',
            subtitle = paste('distance = ', metric, '; transformation = ', abund, sep = ''))
  
  ggsave(file.path(path_fig, paste(metric, abund, 'ind_clone_relatedness.png', sep = '_')), last_plot(), height = 5, width = 8)
  ggsave(file.path(path_fig, paste(metric, abund, 'ind_clone_relatedness.pdf', sep = '_')), last_plot(), height = 5, width = 8)
  
  # run an Adonis test and betadisper
  mod_div_1 <- vegan::adonis(ps_dist ~ evol_fac, data = d_samp, n_perm = 9999)
  mod_betadisper_1 <- betadisper(ps_dist, d_samp$evol_fac)
  
  # save this output
  cat(paste('\n#############################\n# Analysis of individual clones, difference between pre-adaptation history\n#############################\n', sep = ' '), file = file, append = TRUE)
  cat(paste('\n# permutational ANOVA\n', sep = ' '), file = file, append = TRUE)
  print(mod_div_1)
  cat(paste('\n# Homogeneity of variances\n', sep = ' '), file = file, append = TRUE)
  print(anova(mod_betadisper_1))
  
  # this is for what we have relative fitness assays for
  
  # grab the data
  # grab centroids and other data
  d_1 <- get_betadisper_data(mod_betadisper_1)
  
  # combine centroid and eigenvector dataframes for plotting
  betadisper_lines1 <- merge(select(d_1$centroids, group, PCoA1, PCoA2), select(d_1$eigenvector, group, PCoA1, PCoA2), by = c('group'))
  
  # add distances to eigenvector and lines data
  betadisper_lines1 <- mutate(betadisper_lines1, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
  d_1$eigenvector$distances <- d_1$distances$distances
  
  d_1_PCs <- merge(d_1$eigenvector, tibble::rownames_to_column(d_samp, 'sample') %>% select(., sample, pop2, density_ml, fitness), by = 'sample')
  
  p1 <- gather(d_1_PCs, 'trait', 'value', c(density_ml, fitness)) %>%
    gather(., 'PC', 'value2', starts_with('PCoA')) %>%
    filter(., PC %in% c('PCoA1', 'PCoA2')) %>%
    ggplot(., aes(value2, value, col = group)) +
    facet_wrap(~PC + trait, scales = 'free', ncol = 2, labeller = labeller(.multi_line = FALSE)) +
    geom_text(aes(label = pop2)) +
    theme_bw(base_size = 16) +
    scale_color_manual(values = c('darkgrey', 'black')) +
    ggtitle('Correlation between PCoA and measured traits',
            subtitle = paste('distance = ', metric, '; transformation = ', abund, sep = ''))
  
  # save plot, other ways are available
  ggsave(file.path(path_fig, paste(metric, abund, 'regress_traits_vs_PCoA.png', sep = '_')), p1, height = 10, width = 12)
  ggsave(file.path(path_fig, paste(metric, abund, 'regress_traits_vs_PCoA.pdf', sep = '_')), p1, height = 10, width = 12)
  
  # 2. 4 related clones ####
  
  cat(paste('\n#############################\n# Analysis of four clone mixes, difference between pre-adaptation history\n#############################\n', sep = ' '), file = file, append = TRUE)
  
  # filter treatments to keep
  to_keep <- filter(meta_new, treatment %in% c('4_related_clones'))
  ps_sub <- prune_samples(to_keep$SampleID, ps_transform)
  
  # make a data frame of the sample data
  d_samp <- data.frame(sample_data(ps_sub))
  d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                   nclones_fac = as.factor(nclones_fac),
                   evol_fac = as.factor(evolution)) %>%
    column_to_rownames(., 'SampleID')
  
  # calculate distance matrix
  ps_dist <- phyloseq::distance(ps_sub, method = metric)
  
  # run an Adonis test
  mod_div_4 <- vegan::adonis(ps_dist ~ evol_fac, data = d_samp, n_perm = 9999)
  mod_betadisper_4 <- betadisper(ps_dist, d_samp$evol_fac)
  
  # save output
  cat(paste('\n# permutational ANOVA\n', sep = ' '), file = file, append = TRUE)
  print(mod_div_4)
  cat(paste('\n# Homogeneity of variances\n', sep = ' '), file = file, append = TRUE)
  print(anova(mod_betadisper_4))
  
  
  # 3. 24 related clones ####
  
  cat(paste('\n#############################\n# Analysis of 24 clone mixes, difference between pre-adaptation history\n#############################\n', sep = ' '), file = file, append = TRUE)
  
  # filter samples to keep
  to_keep <- filter(meta_new, treatment %in% c('evolved_with_community', 'evolved_without_community'))
  ps_sub <- prune_samples(to_keep$SampleID, ps_transform)
  
  # do an adonis
  # make a data frame of the sample data
  d_samp <- data.frame(sample_data(ps_sub))
  d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                   nclones_fac = as.factor(nclones_fac),
                   evol_fac = as.factor(evolution)) %>%
    column_to_rownames(., 'SampleID')
  
  # calculate distance matrix
  ps_dist <- phyloseq::distance(ps_sub, method = metric)
  
  # run an Adonis test
  mod_div_24 <- vegan::adonis(ps_dist ~ evol_fac, data = d_samp, n_perm = 9999)
  mod_betadisper_24 <- betadisper(ps_dist, d_samp$evol_fac)
  
  # save output
  cat(paste('\n# permutational ANOVA\n', sep = ' '), file = file, append = TRUE)
  print(mod_div_24)
  cat(paste('\n# Homogeneity of variances\n', sep = ' '), file = file, append = TRUE)
  print(anova(mod_betadisper_24))
  
  # Figure: Looking at effect of pre-adaptation history across different levels of diversity ####
  
  # make one big plot of all clones
  to_keep <- filter(meta_new, ! treatment %in% c('4_unrelated_clones'))
  ps_sub <- prune_samples(to_keep$SampleID, ps_transform)
  
  # transform counts to relative abundances for ordination
  ps_sub <- transform_sample_counts(ps_sub, function(x){x / sum(x)})
  
  # wrangle the metadata
  d_samp <- data.frame(sample_data(ps_sub))
  d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                   nclones_fac = as.factor(nclones_fac),
                   evol_fac = as.factor(evolution)) %>%
    column_to_rownames(., 'SampleID') %>%
    unite(., 'id', nclones_fac, evol_fac, sep = ':')
  
  # calculate distance matrix
  ps_dist <- phyloseq::distance(ps_sub, method = metric)
  
  # run a betadisper
  mod_betadisper <- betadisper(ps_dist, d_samp$id)
  
  # grab centroids and other data
  d_fig_preadapt <- get_betadisper_data(mod_betadisper)
  
  # combine centroid and eigenvector dataframes for plotting
  betadisper_lines <- merge(select(d_fig_preadapt$centroids, group, PCoA1, PCoA2), select(d_fig_preadapt$eigenvector, group, PCoA1, PCoA2), by = c('group'))
  
  # add distances to eigenvector and lines data
  betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
  d_fig_preadapt$eigenvector$distances <- d_fig_preadapt$distances$distances
  
  d_fig_preadapt$eigenvalue <- mutate(d_fig_preadapt$eigenvalue, perc = round(eig/sum(eig) * 100, 1))
  
  # split up group into clones and evolution context
  betadisper_lines <- separate(betadisper_lines, group, c('nclones', 'evol'), sep =':') %>%
    mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                               nclones == 'C_high' ~ 'negative_control',
                               TRUE ~ evol),
           nclones =forcats::fct_relevel(nclones,
                                         "C_1", "C_4", "C_24"))
  d_fig_preadapt$centroids <- separate(d_fig_preadapt$centroids, group, c('nclones', 'evol'), sep =':') %>%
    mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                               nclones == 'C_high' ~ 'negative_control',
                               TRUE ~ evol),
           nclones =forcats::fct_relevel(nclones,
                                         "C_1", "C_4", "C_24"))
  d_fig_preadapt$eigenvector <- separate(d_fig_preadapt$eigenvector, group, c('nclones', 'evol'), sep =':') %>%
    mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                               nclones == 'C_high' ~ 'negative_control',
                               TRUE ~ evol),
           nclones =forcats::fct_relevel(nclones,
                                         "C_1", "C_4", "C_24"))
  
  # clone labels
  facet <- c(C_1 = '(a) single clone', C_4 = '(b) 4 clones', C_24 = '(c) 24 clones')
  
  # plot PCoA
  fig_preadapt <- ggplot() +
    # add negative control and lacz ancestor into background
    geom_point(aes(PCoA1, PCoA2, alpha = 1 - distances), select(filter(d_fig_preadapt$eigenvector, evol %in% c('negative_control')), -nclones), size = 0.75, col = 'blue', alpha = 0.2) +
    geom_point(aes(PCoA1, PCoA2), select(filter(d_fig_preadapt$centroids, evol %in% c('negative_control')), -nclones), size = 5, col = 'blue', alpha = 0.2) +
    geom_point(aes(PCoA1, PCoA2, alpha = 1 - distances), select(filter(d_fig_preadapt$eigenvector, evol %in% c('lacz_ancestor')), -nclones), size = 0.75, col = 'orange', alpha = 0.2) +
    geom_point(aes(PCoA1, PCoA2), select(filter(d_fig_preadapt$centroids, evol %in% c('lacz_ancestor')), -nclones), size = 5, col = 'orange', alpha = 0.2) +
    # add points from preadapted clone treatments
    geom_point(aes(PCoA1, PCoA2, col = evol, alpha = 1 - distances), filter(d_fig_preadapt$eigenvector, ! evol %in% c('lacz_ancestor', 'negative_control')), size = 1) +
    geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(filter(betadisper_lines, ! evol %in% c('lacz_ancestor', 'negative_control'))), col = evol, alpha = 1 - distances), filter(betadisper_lines, ! evol %in% c('lacz_ancestor', 'negative_control'))) +
    geom_point(aes(PCoA1, PCoA2, col = evol), filter(d_fig_preadapt$centroids, ! evol %in% c('lacz_ancestor', 'negative_control')), size = 7) +
    theme_bw(base_size = 14, base_family = 'Helvetica') +
    ylab(paste('PCoA Axis 2 (', d_fig_preadapt$eigenvalue$perc[2], ' %)', sep = '')) +
    xlab(paste('PCoA Axis 1 (', d_fig_preadapt$eigenvalue$perc[1], ' %)', sep = '')) +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0)) +
    facet_wrap(~ nclones, labeller = labeller(nclones = facet)) +
    scale_alpha(range = c(0.0001, 1), guide = FALSE) +
    scale_color_manual('', values = c('dark grey', 'black')) +
    ggtitle('Effect of pre-adaptation history on effect of focal species on community composition',
            subtitle = paste('distance = ', metric, '; transformation = ', abund, sep = ''))
    
  
  # save plot, other ways are available
  ggsave(file.path(path_fig, paste(metric, abund, 'effect_of_preadapt_history.png', sep = '_')), fig_preadapt, height = 5, width = 12)
  ggsave(file.path(path_fig, paste(metric, abund,'effect_of_preadapt_history.pdf', sep = '_')), fig_preadapt, height = 5, width = 12)
  
  #-----------------------------------------------------------------------#
  # distance metric just on diversity, ignoring preadaptation history  ####
  #-----------------------------------------------------------------------#
  
  cat(paste('\n#############################\n# Analysis of all clones mixes, ignoring pre-adaptation history\n#############################\n', sep = ' '), file = file, append = TRUE)
  
  # change some values of nclones
  meta_new <- sample_data(ps_transform) %>% data.frame() %>%
    mutate(., n_clones = paste('C_', n_clones, sep = '')) %>%
    mutate(., n_clones = ifelse(treatment == 'lacz_ancestor', 'lacz_ancest', n_clones))
  row.names(meta_new) <- meta_new$SampleID
  sample_data(ps_transform) <- sample_data(meta_new)
  
  # filter samples
  to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0',  'wt_ancestor', 'negative_control'))
  ps_sub <- prune_samples(to_keep$SampleID, ps_transform)
  
  # get the distance matrix out of the data
  ps_dist <- phyloseq::distance(ps_sub, method = metric)
  
  # make a data frame of the sample data
  d_samp <- data.frame(sample_data(ps_sub))
  d_samp <- mutate(d_samp,
                   nclones_fac = as.factor(n_clones)) %>%
    column_to_rownames(., 'SampleID')
  
  # run an Adonis test
  
  mod_nclonefac <-  vegan::adonis(ps_dist ~ nclones_fac, data = d_samp, n_perm = 9999)
  
  # run a multiple comparison to see which treatments are different
  mult_comp <- mctoolsr::calc_pairwise_permanovas(ps_dist, d_samp, 'nclones_fac', n_perm = 9999) %>%
    mutate(., pvalHolm = p.adjust(pval, method = 'holm'),
           pvalHochberg = p.adjust(pval, method = 'hochberg'),
           pvalHommel = p.adjust(pval, method = 'hommel'))
  
  mult_comp_to_save <- mutate(mult_comp, X1 = forcats::fct_recode(X1, `single clone` = 'C_1',
                                                                  `4 clones` = 'C_4',
                                                                  `24 clones` = 'C_24'),
                              X2 = forcats::fct_recode(X2, `LacZ ancestor` = 'lacz_ancest',
                                                       `4 clones` = 'C_4',
                                                       `24 clones` = 'C_24'),
                              contrast = paste(X1, 'vs.', X2, sep = ' ')) %>%
    select(., contrast, everything(),-c(X1, X2)) %>%
    mutate_at(., vars(2:ncol(.)), function(x) signif(x, 2))
  
  # save multiple comparisons out
  write.csv(mult_comp_to_save, file.path('results/mult_comparisons', paste(metric, abund, 'mult_comp.csv', sep = '_')), row.names = FALSE)
  
  # loads of comparisons. Significant differences will be determined by p value correction.
  
  # overwrite metadata to allow plotting of nclones_fac
  sample_data(ps_sub) <- sample_data(d_samp)
  
  # beta-diversity analysis - look at homogeneity of variances
  mod1_dispers <- betadisper(ps_dist, d_samp$nclones_fac)
  
  # save output
  cat(paste('\n# permutational ANOVA\n', sep = ' '), file = file, append = TRUE)
  print(mod_nclonefac)
  cat(paste('\n# Homogeneity of variances\n', sep = ' '), file = file, append = TRUE)
  print(anova(mod1_dispers))
  cat(paste('\n# Multiple comparisons from the pairwise PERMANOVA\n', sep = ' '), file = file, append = TRUE)
  print(mult_comp_to_save)
  
  # Permutation test for F
  pmod <- permutest(mod1_dispers, pairwise = TRUE)
  
  # Tukey's Honest Significant Differences
  T_HSD <- TukeyHSD(mod1_dispers)
  
  # plot distances to centroid - try and play with alpha
  
  # get betadisper data ####
  betadisper_dat <- get_betadisper_data(mod1_dispers)
  
  # do some transformations on the data
  betadisper_dat$eigenvalue <- mutate(betadisper_dat$eigenvalue, percent = round(eig/sum(eig)*100, 1))
  
  # combine centroid and eigenvector dataframes for plotting
  betadisper_lines <- merge(select(betadisper_dat$centroids, group, PCoA1, PCoA2), select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group'))
  
  # add distances to eigenvector and lines data
  betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
  betadisper_dat$eigenvector$distances <- betadisper_dat$distances$distances
  
  # plot PCoA
  p1 <- ggplot() +
    geom_point(aes(PCoA1, PCoA2, col = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24'))), betadisper_dat$eigenvector, size = 1.5) +
    ggConvexHull::geom_convexhull(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$eigenvector, alpha = 0) +
    geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 7) +
    theme_bw(base_size = 16, base_family = 'Helvetica') +
    ylab(paste('PCoA Axis 2 (', betadisper_dat$eigenvalue$percent[2], ' %)', sep = '')) +
    xlab(paste('PCoA Axis 1 (', betadisper_dat$eigenvalue$percent[1], ' %)', sep = '')) +
    theme(legend.position = 'none') +
    scale_color_viridis_d() +
    ggtitle('Effect of diversity',
          subtitle = paste('distance = ', metric, '; transformation = ', abund, sep = ''))
  
  # plot axis 1
  p2 <- ggplot(betadisper_dat$eigenvector, aes(forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24')), PCoA1, col = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24')), fill = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24')))) +
    geom_hline(aes(yintercept = 0)) +
    MicrobioUoE::geom_pretty_boxplot() +
    geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
    theme_bw(base_size = 16, base_family = 'Helvetica') +
    ylab(paste('PCoA Axis 1 (', betadisper_dat$eigenvalue$percent[1], ' %)', sep = '')) +
    xlab('') +
    theme(legend.position = 'none') +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones'))
  
  # screeplot
  p_scree <- ggplot() +
    geom_line(aes(readr::parse_number(PCoA), eig), betadisper_dat$eigenvalue) +
    geom_point(aes(readr::parse_number(PCoA), eig), shape = 21, fill = 'white', betadisper_dat$eigenvalue, size = 2) +
    theme_bw() +
    ylab('Eigenvale') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab('Principal component') +
    ggtitle('Screeplot from PCoA on diversity effect',
            subtitle = paste('distance = ', metric, '; transformation = ', abund, sep = ''))
  
  p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))
  
  paste(metric, abund,'effect_of_preadapt_history.pdf', sep = '_')
  
  ggsave(file.path(path_fig, paste(metric, abund,'PCoA_plot_diversity.png', sep = '_')), p3, height = 5, width = 13)
  ggsave(file.path(path_fig, paste(metric, abund,'PCoA_plot_diversity.pdf', sep = '_')), p3, height = 5, width = 13)
  ggsave(file.path(path_fig, paste(metric, abund,'scree_plot.png', sep = '_')), p_scree, height = 5, width = 6)
  ggsave(file.path(path_fig, paste(metric, abund,'scree_plot.pdf', sep = '_')), p_scree, height = 5, width = 6)
  
  #---------------------------------------#
  # analysis of allopatry vs. sympatry ####
  #---------------------------------------#
  
  cat(paste('\n#############################\n# Analysis of allopatry vs sympatric populations in 4 clone populations\n#############################\n', sep = ' '), file = file, append = TRUE)
  
  # filter treatments to keep
  to_keep <- filter(meta_new, treatment %in% c('4_related_clones', '4_unrelated_clones'))
  ps_sub <- prune_samples(to_keep$SampleID, ps_transform)

  # make a data frame of the sample data
  d_samp <- data.frame(sample_data(ps_sub))
  d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                   nclones_fac = as.factor(nclones_fac),
                   evol_fac = as.factor(evolution),
                   treatment_fac = as.factor(treatment)) %>%
    column_to_rownames(., 'SampleID')
  
  # calculate distance matrix
  ps_dist <- phyloseq::distance(ps_sub, method = metric)
  
  # run an Adonis test
  mod_allopatry <- vegan::adonis(ps_dist ~ treatment_fac, data = d_samp, n_perm = 9999)
  mod_betadisper <- betadisper(ps_dist, d_samp$treatment_fac)
  
  # save output
  cat(paste('\n# permutational ANOVA\n', sep = ' '), file = file, append = TRUE)
  print(mod_allopatry)
  cat(paste('\n# Homogeneity of variances\n', sep = ' '), file = file, append = TRUE)
  print(anova(mod_betadisper))
  closeAllConnections()
}

# lets try run the bloody function
f_analysis(ps2)

1+1
# set up data frame
d_types <- as_tibble(expand.grid(metric = c('unifrac', 'wunifrac', 'jsd', 'jaccard', 'bray'), abund = c('raw', 'hellinger', 'prop', 'log', 'sqrt'), stringsAsFactors = FALSE))

# run each combination
for(i in 1:nrow(d_types)){
  f_analysis(ps2, metric = d_types$metric[i], abund = d_types$abund[i])
  }

