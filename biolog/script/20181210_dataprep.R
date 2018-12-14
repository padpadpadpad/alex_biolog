# prepping biolog files

rm(list = ls())

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggridges)
library(lme4)
library(patchwork)

# figure path
path_fig <- 'biolog/figs'

# extra functions
read_plate <- function(x, type = c('raw_data', 'meta_data')){
  if(missing(type)){type <- 'raw_data'}
  temp <- readxl::read_excel(x) %>%
    janitor::clean_names() %>%
    dplyr::select(1:13) %>%
    tidyr::gather(., 'well', 'od', 2:13) %>%
    dplyr::mutate(., well = readr::parse_number(well),
                  well = paste(x_1, well, sep = '_'),
                  file = basename(tools::file_path_sans_ext(x))) %>%
    dplyr::select(file, well, od)
  if(type == 'raw_data'){temp <- dplyr::mutate(temp, time = file.mtime(x))}
  
  if(type == 'meta_data'){temp <- rename(temp, treatment = od)}
  return(temp)
}

# source extra functions
source('biolog/script/functions.R')

# load in data ####

# list all files
files <- list.files('biolog/data/20181203', full.names = TRUE)

d <- MicrobioUoE::bind_biolog_all('biolog/data/20181203/T6_8.xlsx', sheets = 'Sheet1')

# read in all files
d <- purrr::map_df(files, MicrobioUoE::bind_biolog_all, sheets = 'Sheet1')

# read in metadata for substrates
meta <- read_plate('biolog/data/biolog_ecoplate_metadata.xlsx', type = 'meta_data') %>%
  select(., -file) %>%
  rename(substrate=treatment)
d <- merge(d, meta, by = 'well')

# read in metadata for treatments
meta2 <- read.csv('biolog/data/20181203_metadata.csv', stringsAsFactors = FALSE) %>%
  gather(., set, sample, starts_with('set')) %>%
  merge(., read.csv('biolog/data/20181203_metadata_treatment.csv', stringsAsFactors = FALSE), by = 'sample', all.x = TRUE) %>%
  mutate(., evolved = case_when(sample == 'M9' ~ 'control',
                                grepl('anc', sample) ~ 'ancestor',
                                TRUE ~ evolved),
         population = ifelse(is.na(population), evolved, population))

# bind all pieces together
d <- separate(d, file, c('tp', 'plate'), sep = '_') %>%
  mutate(., set = case_when(readr::parse_number(well) %in% 1:4 ~ 'set_1',
                            readr::parse_number(well) %in% 5:8 ~ 'set_2',
                            readr::parse_number(well) %in% 9:12 ~ 'set_3'))

d <- merge(d, meta2, by = c('set', 'plate'))

# voila
# save this out
write.csv(d, 'biolog/data/20181203_processed.csv', row.names = FALSE)

# data and metadata together
d <- merge(d, meta, by = 'id')

# load in ancestral data
d_ancest <- MicrobioUoE::bind_biolog_all('biolog/data/20170124_Ancestors_gn2biolog.xlsx', sheets = 'Sheet1') %>%
  mutate(id = id + 50,
         pop = 13)
meta_ancest <- data.frame(id = 51:58, treatment = 'wild_type', stringsAsFactors = FALSE)
d_ancest <- merge(d_ancest, meta_ancest, by = 'id') %>%
  filter(., id > 54)
d <- bind_rows(d, d_ancest) %>%
  filter(., id != 49)

# which columns are substrates
Carb_cols <- colnames(d)[grepl('X', colnames(d))]

# stack
d_stack <- gather_(d, 'C_source', 'OD', Carb_cols) %>%
  mutate(., C_source = as.numeric(gsub('X', '', C_source)))

# filter out blank
blank <- filter(d_stack, id == 50) %>%
  rename(., blank = OD) %>%
  select(., blank, C_source)
d_stack <- filter(d_stack, id != 50) %>%
  merge(., blank, by = 'C_source') %>%
  mutate(., OD_cor = OD - blank)
d_noWT <- filter(d_stack, treatment != 'wild_type')
WT <- filter(d_stack, treatment == 'wild_type') %>%
  rename(., WT = OD_cor) %>%
  select(., C_source, WT)

# add column for ranked mean OD_cor
d_noWT <- group_by(d_noWT, C_source) %>%
  mutate(mean_OD = mean(OD_cor)) %>%
  ungroup() %>%
  merge(., WT, by = 'C_source') %>%
  mutate(., rank = dense_rank(desc(mean_OD)))
d_stack2 <- group_by(d_stack, C_source) %>%
  mutate(mean_OD = mean(OD_cor)) %>%
  ungroup() %>%
  merge(., WT, by = 'C_source') %>%
  mutate(., rank = dense_rank(desc(mean_OD)))


# plot performance across wells, ranked by best performance
plot1a <- group_by(d_stack, id) %>%
arrange(., desc(OD_cor)) %>%
  mutate(., rank = 1:96) %>%
  ggplot(.) +
  geom_line(aes(rank, OD_cor, group = id, col = treatment), alpha = 0.25) +
  stat_summary(aes(rank, OD_cor, col = treatment, group = treatment), fun.y = mean, geom = 'line') +
  scale_color_viridis(discrete = TRUE) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = c(0.9, 0.8)) +
  ylab('optical density') +
  xlab('substrate rank') +
  ggtitle('Substrate rank across populations') +
  guides(col = guide_legend(override.aes = list(alpha = 1)))

# plot performance across well, without ranking by best performance
plot1b <- ggplot(d_stack) +
  geom_line(aes(C_source, OD_cor, group = id, col = treatment), alpha = 0.25) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = 'none') +
  ylab('optical density') +
  xlab('substrate') +
  ggtitle('Performance across substrates') +
  scale_color_viridis(discrete = TRUE) 

plot1 <- gridExtra::grid.arrange(plot1a, plot1b, ncol = 1)

ggsave(file.path(path_fig, 'performance_plot.pdf'), plot1, height = 10, width = 10)
ggsave(file.path(path_fig, 'performance_plot2.pdf'), plot1a + facet_wrap(~treatment) + scale_color_manual(values = rep('black', 3)), height = 5, width = 12)


# Calculate phenotypic variance ####

# Take all the distance matrices of each pairwise combination of clones in each treatment
# average Euclidean distance
pop_dists_df <- filter(d, treatment %in% c('comm', 'no_comm')) %>%
  group_by(., treatment) %>%
  do(pop_dists(x = .)) %>%
  data.frame() %>%
  rename(., clone_i = row, clone_j = col)

# create average phenotypic diversity per population
# this is similar to a PCA and PCoA and betadisper()?
V_P <- group_by(pop_dists_df, treatment) %>%
  summarise(., V_P = mean(value)) %>%
  data.frame()

# only two points here!!!

# calculate V_G - genotypic variance ####
# variance of OD across each genotype averaged over all the environments
V_G <- group_by(d_stack, treatment, C_source) %>%
  summarise(V_G = var(OD_cor)) %>%
  data.frame()
V_G_pop <- group_by(V_G, treatment) %>%
  summarise(V_G = mean(V_G)) %>%
  data.frame()

# calculate V_E - environmental variance
# average variance of each clone across all the environments
V_E <- group_by(d_stack, treatment, id) %>%
  summarise(V_E = var(OD_cor)) %>%
  data.frame()
V_E_pop <- group_by(V_E, treatment) %>%
  summarise(V_E = mean(V_E)) %>%
  data.frame()

# analyses
# Genotypic variance
mod_vg <- lmer(V_G ~ treatment + (1|C_source), V_G)
lsmeans::lsmeans(mod_vg, pairwise ~ treatment)
# NOPE

# Phentypic variance
mod_pg <- lm(V_E ~ treatment, V_E)
lsmeans::lsmeans(mod_pg, pairwise ~ treatment)

# plot genotypic and environmental variance across treatments ####
# plot V_G and V_E ####
V_G_plot <- ggplot(V_G, aes(treatment, V_G)) +
  geom_boxplot(aes(fill = treatment, col = treatment), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(treatment, V_G, col = treatment), shape = 21, fill ='white', position = position_jitter(width = 0.1)) +
  ylab('genotypic variance') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Community', 'No Community', 'Wild Type')) +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle(expression(Genotypic~variance~(V[G]))) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  ylim(c(0, 0.46))

V_E_plot <- ggplot(V_E, aes(treatment, V_E)) +
  geom_boxplot(aes(col = treatment, fill = treatment), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(treatment, V_E, col = treatment), shape = 21, fill ='white', position = position_jitter(width = 0.2)) +
  ylab('environmental variance') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Community', 'No Community', 'Wild Type')) +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle(expression(Environmental~variance~(V[E]))) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)

# plot
plot2 <- gridExtra::grid.arrange(V_G_plot, V_E_plot, ncol = 2)

ggsave(file.path(path_fig, 'geno_enviro_var_plot.pdf'), plot2, height = 5, width = 10)

# plot all carbon sources ####
y_axis_change <- group_by(d_stack2, rank) %>%
  summarise(C_source = unique(C_source)) %>%
  pull(C_source)

ggplot(d_stack2) +
  geom_density_ridges2(aes(x = OD_cor, y = factor(rank), fill = treatment, col = treatment), alpha = 0.5, rel_min_height = 0.01) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_point(aes(x = WT, y = factor(rank)), size = 0.5) +
  theme_bw() +
  scale_y_discrete(labels = y_axis_change) +
  ylab('Carbon source') +
  xlab('Optical Density')

ggsave(file.path(path_fig, 'crazy_ggjoy_plot.pdf'), last_plot(), height = 12, width = 6)

# Calculate G x E interaction for each population ####
# see Barrett et al. 2005 Am Nat and Venail et al. 2008 Nature
# 1. calculate responsiveness - indicates differences in the environmental variances and thus measures diversity of resource exploitation strategies (specialists and generalists)
# sum (V_Gj - V_Gi)^2/(2*n_genotypes(n-genotypes - 1))

# create dataframe for standard deviation per clone across environments
d_sd <- group_by(d_stack, treatment, pop, id) %>%
  summarise(., sd_E = sd(OD)) %>%
  data.frame()

# create 2 copies of this for merging later
sd_j_clone <- rename(d_sd, clone_j = id, sd_j = sd_E) 
sd_i_clone <- rename(d_sd, clone_i = id, sd_i = sd_E)

# create every pairwise combination of 1:n (clones/genotypes) for each population
d_R <- group_by(d_sd, treatment, pop) %>%
  do(data.frame(expand.grid(clone_j = .$id, clone_i = .$id))) %>%
  ungroup() %>%
  filter(., clone_j > clone_i) %>%
  merge(., sd_j_clone, by = c('clone_j', 'treatment', 'pop')) %>%
  merge(., sd_i_clone, by = c('clone_i', 'treatment', 'pop'))

# calculate R for each pairwise combination
d_R <- group_by(d_R, treatment, pop) %>%
  mutate(., R_comb = (sd_j - sd_i)^2/(length(unique(clone_i))*(length(unique(clone_i))-1))) %>%
  ungroup()

# calculate responsiveness for each population
# sum of all the pairwise combinations
d_R_pop <- group_by(d_R, treatment, pop) %>%
  summarise(., R_pop = sum(R_comb)) %>%
  data.frame()

# Plot responsiveness
r_plot <- ggplot(d_R_pop, aes(treatment, R_pop)) +
  geom_point(aes(treatment, R_pop, col = treatment), size = 3) +
  ylab('responsiveness') +
  xlab('Treatment') +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('(a) Responsiveness') +
  scale_colour_viridis(discrete = TRUE)

# not significantly different
# summary(lm(R_pop ~ treatment, d_R_pop))

# calculate inconsistency ####

d <- filter(d, id != 50)

# prep data for calculating correlations
d_pearson <- group_by(d, treatment, pop) %>%
  do(pop_cor(x = ., id = 'id', rows_delete = c('treatment', 'id', 'sheet', 'pop'))) %>%
  data.frame()

# merge dataframe to responsiveness dataframe
d_Inconsist <- merge(d_R, d_pearson, by = c('treatment', 'pop', 'clone_j', 'clone_i')) %>%
  mutate(., i = (sd_j*sd_i*(1-pear_cor))/(length(unique(clone_i))*(length(unique(clone_i))-1)))
d_I_pop <- group_by(d_Inconsist, treatment, pop) %>%
  summarise(., I_pop = sum(i),
            pear_pop = mean(pear_cor)) %>%
  data.frame()

# plot inconsistency
I_plot <- ggplot(d_I_pop, aes(treatment, I_pop)) +
  geom_point(aes(treatment, I_pop, col = treatment), size = 3) +
  ylab('Inconsistency') +
  xlab('Treatment') +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('(b) Inconsistency') +
  scale_color_viridis(discrete = TRUE)

p_V_by_G <- r_plot + I_plot

ggsave(file.path(path_fig, 'responsiveness.pdf'), p_V_by_G, height = 5, width = 10)

summary(lm(I_pop ~ treat, d_I_pop))
 
# try a pca ####
d <- unite(d, 'id2', c(treatment, id), sep = '_', remove = FALSE)

# data set ready for PCA
d2 <- filter(d, id != 50 & id != 49)
d_PCA <- d2 %>%
  select(., starts_with('X'))
row.names(d_PCA) <- d2$id2

# create matrix
Euclid_mat <- dist(d_PCA)

# get variables for PCA
d_vars <- select(d2, id, id2, treatment) %>%
  mutate(., treatment = as.factor(treatment))
row.names(d_vars) <- d2$id2
PCA <- prcomp(d_PCA)
biplot(PCA)

# quick and dirty beta disper model
mod_adonis <- vegan::adonis(d_PCA ~ treatment, d_vars) # yes they have different centroids
mctoolsr::calc_pairwise_permanovas(Euclid_mat, d_vars, 'treatment')

Euclid_mat <- dist(d_PCA)
mod <- vegan::betadisper(Euclid_mat, d_vars$treatment)
anova(mod)
TukeyHSD(mod)

# get betadisper dataframes
betadisper_dat <- get_betadisper_data(mod)

# do some transformations on the data
betadisper_dat$eigenvalue <- mutate(betadisper_dat$eigenvalue, percent = eig/sum(eig))

# add convex hull points
betadisper_dat$chull <- group_by(betadisper_dat$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()

# combine centroid and eigenvector dataframes to plot
betadisper_lines <- merge(select(betadisper_dat$centroids, group, PCoA1, PCoA2), select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group'))

ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 4) +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$eigenvector) +
  geom_path(aes(PCoA1, PCoA2, col = group), betadisper_dat$chull ) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, col = group), betadisper_lines) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('PCoA Axis 2 [10.3%]') +
  xlab('PCoA Axis 1 [36.1%]') +
  scale_color_viridis('', discrete = TRUE) +
  #coord_fixed(sqrt(betadisper_dat$eigenvalue$percent[2]/betadisper_dat$eigenvalue$percent[1])) +
  coord_fixed() +
  theme(legend.position = 'top') +
  ggtitle('PCoA across treatments') +
  guides(col = guide_legend(ncol = 8))

ggsave(file.path(path_fig, 'PCoA_across_treatments.pdf'), last_plot(), height = 5, width = 7)

# distance plot
ggplot(betadisper_dat$distances, aes(group, distances, fill = group, col = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.2)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_viridis('', discrete = TRUE, labels = c('Community', 'No Community', 'Wild Type')) +
  scale_fill_viridis('', discrete = TRUE, labels = c('Community', 'No Community', 'Wild Type')) +
  ylab('Distance to centroid') +
  xlab('') +
  scale_x_discrete(labels = c('Community', 'No Community', 'Wild Type'))

mod <- lm(distances ~ group, betadisper_dat$distances)
summary(mod)   
