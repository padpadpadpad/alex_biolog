# analysis script
# comparing variance in the presence of new ancestor biolog plates

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

# source extra functions
source('biolog/script/functions.R')

# load in data ####
d <- MicrobioUoE::bind_biolog_all('biolog/data/biolog_data_assay_1.xlsx', sheets = 'Sheet1')

# make into long format ####

# meta data
meta <- data.frame(id = 1:50, treatment = c(rep('comm', times = 24), rep('no_comm', times = 24), rep('wild_type', times = 2)), stringsAsFactors = FALSE) %>%
  mutate(., pop = c(rep(1:12, each = 4), 13, 13))

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

# load in repeat biolog plates

# vector of metadata
new_biologs_id <- Hmisc::Cs(ctrl, a10, a9, a8, a7, a6, a5, a4, a3, a2, 32, 29, 17, 11, 6, a1, 48, 44, 42, 35)

d_ancest_2 <- MicrobioUoE::bind_biolog_all('biolog/data/Biolog_gn2-2018.xlsx', sheets = 'Sheet1') %>%
  mutate(., id = new_biologs_id,
         type = 'new_biolog',
         treatment = case_when(id == 'ctrl' ~ 'control',
                               grepl('a', id) ~ 'wild_type'))

# which columns are substrates
Carb_cols <- colnames(d)[grepl('X', colnames(d))]

# compare previous biolog run to the new repeated biolog plates
d_twice <- filter(d, id %in% new_biologs_id) %>%
  mutate(., type = 'old_biolog',
         id = as.character(id)) %>%
  bind_rows(., filter(d_ancest_2, id %in% .$id)) %>%
  group_by(., id) %>%
  mutate(., treatment = ifelse(is.na(treatment), treatment[!is.na(treatment)], treatment),
         pop = ifelse(is.na(pop), pop[!is.na(pop)], pop)) %>%
  ungroup() %>%
  gather(., 'C_source', 'OD', starts_with('X')) %>%
  spread(., type, OD) %>%
  mutate(C_source = readr::parse_number(C_source))

ggplot(d_twice, aes(new_biolog, old_biolog, col = id)) +
  geom_line() +
  ylim(c(0, 3)) +
  xlim(c(0,3)) +
  theme_bw() +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  ylab('OD of old biolog plate') +
  xlab('OD of new biolog plate') +
  ggtitle('Comparison between duplicates IDs of P. fluorescens clones run on both old and new biolog plates',
          subtitle = 'Dashed line would be the expected 1:1 line')

ggsave(file.path(path_fig, 'old_vs_new_biologs.pdf'), last_plot(), height = 5, width = 7)

# the new biologs seem much more productive than the old biologs

# check the comparison between the original blank and the new blank
d_control <- filter(d, id == 50) %>%
  select(., -c(pop)) %>%
  mutate(., type = 'old_biolog',
         id = 'ctrl',
         treatment = 'control') %>%
  bind_rows(., filter(d_ancest_2, id == 'ctrl')) %>%
  gather(., 'C_source', 'OD', starts_with('X'))
  mutate(C_source = readr::parse_number(C_source))

ggplot(d_control, aes(type, OD)) +
  MicrobioUoE::geom_pretty_boxplot(aes(col = type, fill = type)) +
  geom_point(aes(col = type), shape = 21, fill = 'white', position = position_jitter(height = 0, width = 0.1)) +
  theme_bw() +
  xlab('Blank') +
  ggtitle('Difference in blank readings between old and new biolog plates') +
  theme(legend.position = 'none')

ggsave(file.path(path_fig, 'old_vs_new_biologs_control.pdf'), last_plot(), height = 5, width = 6)

# compare the G and E variances including the new ancestors (do not correct for blanks) ####
d_ancest_new <- filter(d_ancest_2, treatment == 'wild_type') %>%
  mutate(., pop = readr::parse_number(id),
         id = 1:10 + max(d$id))

# stack new data and old data
d <- mutate(d, type = 'old_biolog') %>%
  bind_rows(., d_ancest_new)


d_stack <- gather(d, C_source, OD, starts_with('X')) %>%
  mutate(., C_source = readr::parse_number(C_source))

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
arrange(., desc(OD)) %>%
  mutate(., rank = 1:96) %>%
  ggplot(.) +
  geom_line(aes(rank, OD, group = id, col = treatment), alpha = 0.25) +
  stat_summary(aes(rank, OD, col = treatment, group = interaction(treatment, type)), fun.y = mean, geom = 'line') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = c(0.9, 0.8)) +
  ylab('optical density') +
  xlab('substrate rank') +
  ggtitle('Substrate rank across populations') +
  guides(col = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ treatment)

group_by(filter(d_stack, type == 'new_biolog'), id) %>%
  arrange(., desc(OD)) %>%
  mutate(., rank = 1:96) %>%
  ggplot(.) +
  geom_line(aes(C_source, OD, col = id), alpha = 0.25) +
  scale_color_viridis() +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = c(0.9, 0.8)) +
  ylab('optical density') +
  xlab('substrate rank') +
  ggtitle('Substrate rank across populations') +
  guides(col = guide_legend(override.aes = list(alpha = 1)))

# plot performance across well, without ranking by best performance
plot1b <- ggplot(d_stack) +
  geom_line(aes(C_source, OD, group = id, col = treatment), alpha = 0.25) +
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
V_G <- group_by(d_stack, treatment, C_source, type) %>%
  summarise(V_G = var(OD_cor)) %>%
  data.frame()
V_G_pop <- group_by(V_G, treatment, type) %>%
  summarise(V_G = mean(V_G)) %>%
  data.frame()

# calculate V_E - environmental variance
# average variance of each clone across all the environments
V_E <- group_by(d_stack, treatment, id, type) %>%
  summarise(V_E = var(OD_cor)) %>%
  data.frame()
V_E_pop <- group_by(V_E, treatment, type) %>%
  summarise(V_E = mean(V_E)) %>%
  data.frame()

# analyses
# Genotypic variance
mod_vg <- lmer(V_G ~ treatment * type + (1|C_source), V_G)
lsmeans::lsmeans(mod_vg, pairwise ~ treatment*type)
# NOPE

# Phentypic variance
mod_pg <- lm(V_E ~ treatment, V_E)
lsmeans::lsmeans(mod_pg, pairwise ~ treatment)

# plot genotypic and environmental variance across treatments ####
# plot V_G and V_E ####
V_G_plot <- ggplot(V_G, aes(interaction(treatment, type), V_G)) +
  geom_boxplot(aes(fill = treatment, col = treatment), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(interaction(treatment, type), V_G, col = treatment), shape = 21, fill ='white', position = position_jitter(width = 0.1)) +
  ylab('genotypic variance') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Wild Type New', 'Community', 'No Community', 'Wild Type Old')) +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle(expression(Genotypic~variance~(V[G]))) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  ylim(c(0, 0.46))

V_E_plot <- ggplot(V_E, aes(interaction(treatment, type), V_E)) +
  geom_boxplot(aes(col = treatment, fill = treatment), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(interaction(treatment, type), V_E, col = treatment), shape = 21, fill ='white', position = position_jitter(width = 0.2)) +
  ylab('environmental variance') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Wild Type New', 'Community', 'No Community', 'Wild Type Old')) +
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
d <- mutate(d, treat_2 = case_when(treatment == 'wild_type' ~ paste(treatment, type, sep = '_'),
                                   TRUE ~ treatment)) %>%
  unite(., 'id2', c(treat_2, id), sep = '_', remove = FALSE)

# data set ready for PCA
d2 <- filter(d, id != 50 & id != 49)
d_PCA <- d2 %>%
  select(., starts_with('X'))
row.names(d_PCA) <- d2$id2

# create matrix
Euclid_mat <- dist(d_PCA)

# get variables for PCA
d_vars <- select(d2, id, id2, treat_2) %>%
  mutate(., treatment = as.factor(treat_2))
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

p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 4) +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$eigenvector) +
  geom_path(aes(PCoA1, PCoA2, col = group), betadisper_dat$chull ) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, col = group), betadisper_lines) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab('PCoA Axis 2 [0.029%]') +
  xlab('PCoA Axis 1 [85.4%]') +
  scale_color_viridis('', discrete = TRUE) +
  #coord_fixed(sqrt(betadisper_dat$eigenvalue$percent[2]/betadisper_dat$eigenvalue$percent[1])) +
  coord_fixed() +
  theme(legend.position = 'top') +
  ggtitle('PCoA across treatments') +
  guides(col = guide_legend(ncol = 8))

ggsave(file.path(path_fig, 'PCoA_across_treatments.pdf'), last_plot(), height = 5, width = 7)

# distance plot
p2 <- ggplot(betadisper_dat$distances, aes(group, distances, fill = group, col = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.2)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_viridis('', discrete = TRUE, labels = c('Community', 'No Community', 'Wild Type New', 'Wild Type Old')) +
  scale_fill_viridis('', discrete = TRUE, labels = c('Community', 'No Community', 'Wild Type New', 'Wild Type Old')) +
  ylab('Distance to centroid') +
  xlab('') +
  scale_x_discrete(labels = c('Community', 'No Community', 'Wild Type New', 'Wild Type Old')) +
  ggtitle('Distance to centroid of each treatment based on a PCoA')

ggsave(file.path(path_fig, 'dist_2_centroid.pdf'), last_plot(), height = 5, width = 7)


mod <- lm(distances ~ group, betadisper_dat$distances)
summary(mod)   
