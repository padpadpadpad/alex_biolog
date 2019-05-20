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
library(flextable)
library(webshot)
library(rmarkdown)
library(officer)
# if not installed, install mctoolsr run devtools::install_github('leffj/mctoolsr')

#--------------#
# functions ####
#--------------#

# convert distance matrix to dataframe
dist_2_df <- function(dist_ob){
  m <- as.matrix(dist_ob) # coerce dist object to a matrix
  xy <- t(combn(colnames(m), 2))
  return(data.frame(xy, dist=m[xy], stringsAsFactors = FALSE))
}

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
sample_sums(ps2)

# rarefy once
#ps_rare <- rarefy_even_depth(ps2)
#saveRDS(ps_rare, 'sequencing/data/output/20171024_17:18/ps_rarefied.rds')
ps_rare <- readRDS('sequencing/data/output/20171024_17:18/ps_rarefied.rds')
# ps2 <- ps_rare

# remove Pseudomonas fluorescens reads from the data, leave other pseudomonads in
SBW25 = "ACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAA"

# create distance from P. fluorescens table
#d_otu_phylodist <- adephylo::distTips(phy_tree(ps2))
#d_otu_phylodist <- dist_2_df(d_otu_phylodist)
#d_otu_phylodist <- filter(d_otu_phylodist, X1 == SBW25) %>%
  #select(otu = X2, dist) %>%
  #saveRDS('sequencing/data/output/20171024_17:18/otu_dist_from_SBW25.rds')
  
ps2 <- subset_taxa(ps2, ! rownames(tax_table(ps2)) %in% c(SBW25))

# -------------------------------------------------------------------------------#
# Analysis of the preadaptation with and without nmc on community composition ####
#--------------------------------------------------------------------------------#

# perfectly replicated for 1, 4 and 24 clones
# do these separately and then pvalue correct

# 1. individual clones ####

# still unsure whether I need to nest these and how that would help me
# reps are not paired between treatments so not sure restricted shuffling is necessary

# samples to keep - just individual clones
to_keep <- filter(meta_new, treatment %in% c('individual_clone'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
# weighted Unifrac - relative abundances
# Bray-Curtis - absolute abundances
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# function
fac_to_letters <- function(fac){
  
}

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution),
                 preadapt_pop = as.factor(preadapt_pop),
                 pop2 = preadapt_pop) %>%
  column_to_rownames(., 'SampleID')
levels(d_samp$pop2) <- letters[1:length(levels(d_samp$pop2))]

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# side note
# see whether clones within replicates are on average closer to each other than clones in other replicates

# calculate average distance of each clone to its sympatyric populations and allopatric populations across reps
d_wunifrac <- dist_2_df(ps_wunifrac) %>%
  merge(., select(d_samp, X1_evol = evolution, X1_preadapt_pop = preadapt_pop) %>% rownames_to_column(., var = 'X1'), by = 'X1') %>%
  merge(., select(d_samp, X2_evol = evolution, X2_preadapt_pop = preadapt_pop) %>% rownames_to_column(., var = 'X2'), by = 'X2') %>%
  filter(., X1_evol == X2_evol) %>%
  mutate(same_rep = ifelse(X1_preadapt_pop == X2_preadapt_pop, 'Y', 'N')) %>%
  gather(., 'random', 'clone', c(X1, X2)) %>%
  group_by(clone, same_rep, X1_evol) %>%
  summarise(., mean = mean(dist)) %>%
  ungroup()

# check n - is correct
group_by(d_wunifrac, same_rep, X1_evol) %>%
  tally()

labels <- c(with_community = '(a) pre-adapted with nmc', without_community = '(b) pre-adapted without nmc')

ggplot(d_wunifrac, aes(same_rep, mean)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, size = 3, position = position_jitter(width = 0.1)) +
  facet_wrap(~ X1_evol, labeller = labeller(X1_evol = labels)) +
  theme_bw(base_size = 16) +
  ylab('weighted Unifrac distance') +
  xlab('') +
  scale_x_discrete(labels = c('allopatric pair', 'sympatric pair')) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0))

ggsave(file.path(path_fig, 'ind_clone_relatedness.png'), last_plot(), height = 5, width = 8)
ggsave(file.path(path_fig, 'ind_clone_relatednes.pdf'), last_plot(), height = 5, width = 8)

# run an Adonis test and betadisper
mod_div_1 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_1 <- betadisper(ps_wunifrac, d_samp$evol_fac)

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
  ggtitle('Correlation between PCoA and measured traits',
          subtitle = 'PCoA1 ~ 44% & PCoA2 ~ 21%') +
  scale_color_manual(values = c('darkgrey', 'black'))

# save plot, other ways are available
ggsave(file.path(path_fig, 'regress_traits_vs_PCoA.png'), p1, height = 10, width = 12)
ggsave(file.path(path_fig, 'regress_traits_vs_PCoA.pdf'), p1, height = 10, width = 12)

# 2. 4 related clones ####


# filter treatments to keep
to_keep <- filter(meta_new, treatment %in% c('4_related_clones'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution)) %>%
  column_to_rownames(., 'SampleID')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# run an Adonis test
mod_div_4 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_4 <- betadisper(ps_wunifrac, d_samp$evol_fac)

# 3. 24 related clones ####

# filter samples to keep
to_keep <- filter(meta_new, treatment %in% c('evolved_with_community', 'evolved_without_community'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

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

# run an Adonis test
mod_div_24 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_24 <- betadisper(ps_wunifrac, d_samp$evol_fac)

# Figure: Looking at effect of pre-adaptation history across different levels of diversity ####

# make one big plot of all clones
to_keep <- filter(meta_new, ! treatment %in% c('4_unrelated_clones'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# wrangle the metadata
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
d_fig_preadapt <- get_betadisper_data(mod_betadisper)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(d_fig_preadapt$centroids, group, PCoA1, PCoA2), select(d_fig_preadapt$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
d_fig_preadapt$eigenvector$distances <- d_fig_preadapt$distances$distances

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
  geom_point(aes(PCoA1, PCoA2, alpha = 1 - distances), select(filter(d_fig_preadapt$eigenvector, evol %in% c('lacz_ancestor')), -nclones), size = 0.75, col = '#e41a1c', alpha = 0.4) +
  geom_point(aes(PCoA1, PCoA2), select(filter(d_fig_preadapt$centroids, evol %in% c('lacz_ancestor')), -nclones), size = 5, col = '#e41a1c', alpha = 0.4) +
  # add points from preadapted clone treatments
  geom_point(aes(PCoA1, PCoA2, col = evol, alpha = 1 - distances), filter(d_fig_preadapt$eigenvector, ! evol %in% c('lacz_ancestor', 'negative_control')), size = 1) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(filter(betadisper_lines, ! evol %in% c('lacz_ancestor', 'negative_control'))), col = evol, alpha = 1 - distances), filter(betadisper_lines, ! evol %in% c('lacz_ancestor', 'negative_control'))) +
  geom_point(aes(PCoA1, PCoA2, col = evol), filter(d_fig_preadapt$centroids, ! evol %in% c('lacz_ancestor', 'negative_control')), size = 7) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ylab('PCoA Axis 2 (11.7 %)') +
  xlab('PCoA Axis 1 (27.9 %)') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  facet_wrap(~ nclones, labeller = labeller(nclones = facet)) +
  scale_alpha(range = c(0.0001, 1), guide = FALSE) +
  scale_color_manual('', values = c('dark grey', 'black'))

# save plot, other ways are available
ggsave(file.path(path_fig, 'effect_of_preadapt_history.png'), fig_preadapt, height = 5, width = 12)
ggsave(file.path(path_fig, 'effect_of_preadapt_history.pdf'), fig_preadapt, height = 5, width = 12)

#---------------------------------------#
# analysis of allopatry vs. sympatry ####
#---------------------------------------#

# filter treatments to keep
to_keep <- filter(meta_new, treatment %in% c('4_related_clones', '4_unrelated_clones'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution),
                 treatment_fac = as.factor(treatment)) %>%
  column_to_rownames(., 'SampleID')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# run an Adonis test
mod_allopatry <- vegan::adonis(ps_wunifrac ~ treatment_fac, data = d_samp, n_perm = 9999)
mod_betadisper <- betadisper(ps_wunifrac, d_samp$treatment_fac)

plot(mod_betadisper)

# not significant but need to make the plot

# make table of effect of pre-adaptation ####
d_table <- tibble(clone = c(1, 4, 24), mod = list(mod_div_1, mod_div_4, mod_div_24))

tidy_adonis <- function(adonis_mod){
  stuff <- adonis_mod$aov.tab %>% data.frame(stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column(var = 'factor') %>%
    janitor::clean_names() %>%
    rename(p_value = pr_f)
  return(stuff)
}

d_table <- mutate(d_table, output = purrr::map(mod, tidy_adonis)) %>%
  unnest(output) %>%
  filter(., factor != 'Total') %>%
  spread(factor, df) %>%
  fill(., Residuals, .direction = 'up') %>%
  filter(!is.na(p_value)) %>%
  unite(., 'd.f.', c(evol_fac, Residuals), sep = ', ') %>%
  select(clone, f_model, `d.f.`, r2, p_value) %>%
  mutate_at(., vars(f_model, r2), function(x) round(x, 2)) %>%
  mutate_all(., as.character)

# super_fp
super_fp <- fp_text(vertical.align = "superscript", font.size = 8, font.family = 'Times')

# italics_fp
italic_fp <- fp_text(italic = TRUE, font.size = 16, font.family = 'Times')

table <- flextable(d_table) %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(clone = 'Clonal diversity',
                    f_model = "F statistic",
                    p_value = "p value") %>%
  compose(., j = "r2", part = "header", 
          value = as_paragraph("R", as_chunk("2", props = super_fp))) %>%
  compose(., j = "d.f.", part = "header", 
          value = as_paragraph(as_chunk("d.f.", props = italic_fp))) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit()

# save as a png ####
# create an Rmd file
rmd_name <- tempfile(fileext = ".Rmd")
cat("```{r echo=FALSE}\ntable\n```", file = rmd_name)

# render as an html file ----
html_name <- tempfile(fileext = ".html")
render(rmd_name, output_format = "html_document", output_file = html_name )

# get a png from the html file with webshot ----
webshot(html_name, zoom = 2, file = "table_2.png", 
        selector = "body > div.container-fluid.main-container > div.tabwid > table")

#-----------------------------------------------------------------------#
# Weighted Unifrac just on diversity, ignoring preadaptation history ####
#-----------------------------------------------------------------------#

# change some values of nclones
meta_new <- sample_data(ps2) %>% data.frame() %>%
  mutate(., n_clones = paste('C_', readr::parse_number(n_clones), sep = '')) %>%
  mutate(., n_clones = case_when(treatment == 'lacz_ancestor' ~ 'lacz_ancest',
                                 treatment == 'negative_control' ~ 'C_high',
                                 TRUE ~ as.character(n_clones)))
row.names(meta_new) <- meta_new$SampleID
sample_data(ps2) <- sample_data(meta_new)

# filter samples
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0',  'wt_ancestor', 'negative_control'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# make counts proportions
ps_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})
# hellinger proportions
#ps_prop <- transform_sample_counts(ps_sub, function(x){sqrt(x / sum(x))})

# define distance metric
metric = 'wunifrac'

# get the distance matrix out of the data
ps_dist <- phyloseq::distance(ps_prop, method = metric)

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp,
                 nclones_fac = as.factor(n_clones)) %>%
  column_to_rownames(., 'SampleID')

# run an Adonis test
mod_nclonefac <- vegan::adonis(ps_dist ~ nclones_fac, data = d_samp, n_perm = 9999)

# run a multiple comparison to see which treatments are different
mult_comp <- mctoolsr::calc_pairwise_permanovas(ps_dist, d_samp, 'nclones_fac', n_perm = 9999)

d_num <- group_by(d_samp, nclones_fac) %>% tally()

mult_comp_to_save <- merge(mult_comp, rename(d_num, X1 = nclones_fac, n_1 = n), by = 'X1') %>%
  merge(., rename(d_num, X2 = nclones_fac, n_2 = n), by = 'X2') %>%
  unite(., 'n', c(n_1, n_2), sep = ', ') %>%
  mutate(., X1 = forcats::fct_recode(X1, `single clone` = 'C_1',
                                                                `4 clones` = 'C_4',
                                                                `24 clones` = 'C_24',
                                     `negative control` = 'C_high'),
                            X2 = forcats::fct_recode(X2, `LacZ ancestor` = 'lacz_ancest',
                                                     `4 clones` = 'C_4',
                                                     `24 clones` = 'C_24',
                                                     `negative control` = 'C_high'),
                            contrast = paste(X1, 'vs.', X2, sep = ' ')) %>%
  select(., contrast, n, everything(),-c(X1, X2)) %>%
  mutate_at(., vars(3:ncol(.)), function(x) signif(x, 2)) %>%
  arrange(desc(contrast))

# make table
super_fp <- officer::fp_text(vertical.align = "superscript", font.size = 8, font.family = 'Times')
sub_fp <- officer::fp_text(vertical.align = "subscript", font.size = 8, font.family = 'Times')
italic_fp <- officer::fp_text(italic = TRUE, font.size = 16, font.family = 'Times')

table <- flextable(select(mult_comp_to_save, contrast, n, R2, pval, pvalFDR)) %>%
  align(j = c('contrast'), align = 'left', part = 'all') %>%
  align(j = c('R2', 'pval', 'pvalFDR', 'n'), align = 'center', part = 'all') %>%
  add_footer_row(., colwidths = 5, values = '') %>% 
  bold(i = ~ pval < 0.05, j = ~ pval) %>%
  bold(i = ~ pvalFDR < 0.05, j = ~ pvalFDR) %>%
  set_header_labels(n = 'number of samples per treatment',
                    pval = "raw p value") %>%
  compose(., j = "R2", part = "header", 
          value = as_paragraph("R", as_chunk("2", props = super_fp))) %>%
  compose(., j = "pvalFDR", part = "header", 
          value = as_paragraph("p", as_chunk("adj", props = sub_fp))) %>%
  compose(., j = "contrast", part = "footer", 
          value = as_paragraph("p", as_chunk("adj", props = sub_fp), 'was calculated using the false discovery rate method')) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  padding(padding.top = 5, part = 'footer')

# save as a png ####
# create an Rmd file
rmd_name <- tempfile(fileext = ".Rmd")
cat("```{r echo=FALSE}\ntable\n```", file = rmd_name)

# render as an html file ----
html_name <- tempfile(fileext = ".html")
render(rmd_name, output_format = "html_document", output_file = html_name )

# get a png from the html file with webshot ----
webshot(html_name, zoom = 2, file = "table_3.png", 
        selector = "body > div.container-fluid.main-container > div.tabwid > table")

# save multiple comparisons out
write.csv(mult_comp_to_save, paste('sequencing/data/output/',metric, '_mult_comp.csv', sep = ''), row.names = FALSE)

# loads of comparisons. Significant differences will be determined by p value correction.

# overwrite metadata to allow plotting of nclones_fac
sample_data(ps_prop) <- sample_data(d_samp)

# beta-diversity analysis - look at homogeneity of variances
mod1_dispers <- betadisper(ps_dist, d_samp$nclones_fac)

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
  geom_point(aes(PCoA1, PCoA2, col = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), alpha = 0.5 - distances), betadisper_dat$eigenvector, size = 1.5) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), alpha = 0.5 - distances), betadisper_lines) +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 7) +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab('PCoA Axis 2 (18.9%)') +
  xlab('PCoA Axis 1 (44.6%)') +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  ggtitle('(a)')

# plot axis 1
p2 <- ggplot(betadisper_dat$eigenvector, aes(forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), PCoA1, col = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), fill = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')))) +
  geom_hline(aes(yintercept = 0)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab('PCoA Axis 1 (44.6%)') +
  xlab('') +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones', 'negative\ncontrol')) +
  ggtitle('(b)')

group_by(betadisper_dat$eigenvector, group) %>%
  summarise(above_1 = sum(PCoA1 > 0)/n())
  

# screeplot
p_scree <- ggplot() +
  geom_line(aes(readr::parse_number(PCoA), eig), betadisper_dat$eigenvalue) +
  geom_point(aes(readr::parse_number(PCoA), eig), shape = 21, fill = 'white', betadisper_dat$eigenvalue, size = 2) +
  theme_bw() +
  ylab('Eigenvale') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab('Principal component')

p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))

ggsave(file.path(path_fig, 'PCoA_plot_diversity.png'), p3, height = 6, width = 12)
ggsave(file.path(path_fig, 'PCoA_plot_diversity.pdf'), p3, height = 6, width = 12)
ggsave(file.path(path_fig, 'scree_plot.png'), p_scree, height = 5, width = 6)
ggsave(file.path(path_fig, 'scree_plot.pdf'), p_scree, height = 5, width = 6)

# what correlates with PC1?

# first approach - calculate abundance change with PCoA1 with each sample

d_otu_phylodist <- readRDS('sequencing/data/output/20171024_17:18/otu_dist_from_SBW25.rds')

d_ps <- psmelt(ps_prop) %>%
  merge(., select(betadisper_dat$eigenvector, Sample = sample, pcoa1 = PCoA1), by = 'Sample') %>%
  janitor::clean_names() %>%
  merge(., d_otu_phylodist, by = 'otu')

d_ps2 <- group_by(d_ps, kingdom, phylum, class, sample, n_clones, pcoa1) %>%
  summarise(., abundance = sum(abundance),
            dist = mean(dist)) %>%
  ungroup() %>%
  filter(!is.na(class)) %>%
  nest(-class) %>%
  mutate(., model = purrr::map(data, ~cor.test(~ abundance + pcoa1, .x) %>% broom::tidy())) %>%
  unnest(model) %>%
  mutate(., padj = p.adjust(p.value, method = 'fdr')) %>%
  filter(., padj < 0.05 & !is.na(padj)) %>%
  unnest(data)

ggplot(d_ps2, aes(pcoa1, abundance, col = dist)) +
  geom_point() +
  stat_smooth(aes(group = class), method = 'lm', se = FALSE) +
  facet_wrap(~ phylum)

ggplot(d_ps2, aes(dist, estimate)) +
  geom_point()

d_ps <- tax_glom(ps_prop, taxrank = 'Phylum', NArm = TRUE) %>% 
  psmelt(.) %>%
  merge(., select(betadisper_dat$eigenvector, Sample = sample, pcoa1 = PCoA1), by = 'Sample') %>%
  janitor::clean_names()


# agglomerate taxa at the Genus level
ps_prop2 <- tax_glom(ps_prop, taxrank = 'Genus', NArm = TRUE)
ps_sub2 <- tax_glom(ps_sub, taxrank = 'Genus', NArm = TRUE)
ps_dist <- distance(ps_prop2, 'bray')

nmds <- vegan::metaMDS(ps_dist, distance = 'bray')
stressplot(nmds)
plot(nmds)

ps_counts <- otu_table(ps_sub2) %>%
  data.frame()

# do a double principle component analysis
ps_dpcoa <- DPCoA(ps_prop2)

d_taxa <- psmelt(ps_prop2) %>%
  janitor::clean_names() %>%
  select(otu, kingdom, phylum, class, order, family, genus) %>%
  distinct_all()

# grab data out 
d_dpcoa_sample <- tibble(axis_1_scores = ps_dpcoa$li$Axis1,
                         axis_2_scores = ps_dpcoa$li$Axis2,
                         type = 'sample',
                         sample = row.names(ps_dpcoa$li))
d_dpcoa_species <- tibble(axis_1_scores = ps_dpcoa$dls$CS1, 
                          axis_2_scores = ps_dpcoa$dls$CS2,
                          type = 'species',
                          otu = row.names(ps_dpcoa$dls)) %>%
  merge(., d_taxa, by = 'otu', all.x = TRUE)

d_dpcoa_sample <- merge(d_dpcoa_sample, rownames_to_column(d_samp, 'sample') %>% select(sample, n_clones), by = 'sample')

p1 <- gather(d_dpcoa_sample, 'axis', 'score', starts_with('axis')) %>%
  ggplot(., aes(forcats::fct_relevel(n_clones, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), score, col = forcats::fct_relevel(n_clones, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), fill = forcats::fct_relevel(n_clones, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')))) +
  geom_hline(aes(yintercept = 0)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab('DPCoA Axis Score') +
  xlab('') +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones', 'negative\ncontrol')) +
  facet_wrap(~axis, labeller = labeller(axis = c(axis_1_scores = 'Axis 1: 48.6 %', axis_2_scores = 'Axis 2: 25.3 %')))

p2 <- gather(d_dpcoa_species, 'axis', 'score', starts_with('axis')) %>%
  ggplot(., aes(forcats::fct_reorder(phylum, score, median), score)) +
  geom_hline(aes(yintercept = 0)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab('DPCoA Axis Score') +
  xlab('Phylum') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~axis, labeller = labeller(axis = c(axis_1_scores = 'Axis 1: 48.6 %', axis_2_scores = 'Axis 2: 25.3 %')))

p1 + p2 + plot_layout(ncol = 1)

# grab 

p1 <- plot_ordination(ps_prop2, ps_dpcoa, type = 'samples', col = 'n_clones') +
  stat_ellipse()
p2 <- plot_ordination(ps_prop2, ps_dpcoa, type = 'species', col = 'Phylum')

p1 + p2

d_ps2 <- filter(d_ps, Abundance > 0) %>%
  group_by(OTU) %>%
  filter(., n() > 10) %>%
  nest(-OTU) %>%
  mutate(., model = purrr::map(data, ~cor.test(~ Abundance + PCoA1, .x) %>% broom::tidy())) %>%
  unnest(model) %>%
  filter(., p.value < 0.05) %>%
  unnest(data)
  


ggplot(filter(d_ps2, Genus == 'Chitinophaga'), aes(PCoA1, Abundance)) +
  geom_point() +
  facet_wrap(~OTU, scales = 'free_y')
