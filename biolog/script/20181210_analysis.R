# analysis script

rm(list = ls())

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggridges)
library(lme4)
library(patchwork)
library(widyr)
library(corrr)

# figure path
path_fig <- 'biolog/figs'

# source extra functions
source('biolog/script/functions.R')

# load in data ####
# read in all files
d <- read.csv('biolog/data/20181203_processed.csv', stringsAsFactors = FALSE)

# look at variation across clones

#------------------------------------------------#
# check blanks and create an OD cor dataframe ####
#------------------------------------------------#

# filter out blank
blank <- filter(d, sample == 'M9')

ggplot(blank, aes(substrate, od, col = as.factor(od_wave))) +
  geom_point() +
  facet_wrap(~tp) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

blank_ave <- group_by(blank, substrate, od_wave) %>%
  summarise(., blank = median(od),
            sd_blank = sd(od))

d <- filter(d, sample != 'M9') %>%
  merge(., blank_ave, by = c('substrate', 'od_wave')) %>%
  mutate(., od_cor = od - blank)

# look at variation across clone reps
d_sameclone <- filter(d, sample == 'anc1')
ggplot(d_sameclone, aes(substrate, od_cor, col = interaction(set, plate), shape = as.factor(od_wave))) +
  geom_point() +
  facet_wrap(~tp + od_wave, labeller = labeller(.multi_line = FALSE)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

d_ancest <- filter(d, evolved == 'ancestor')
ggplot(d_ancest, aes(substrate, od_cor, col = sample, shape = as.factor(od_wave))) +
  geom_point(position = position_dodge(width = 1)) +
  facet_wrap(~ tp + od_wave, labeller = labeller(.multi_line = FALSE)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# lets use T4 590
d_t4_590 <- filter(d, od_wave == 590 & tp == 'T4')

# remove anc1 because it is super strange - this is a bit naughty
d_t4_590 <- filter(d_t4_590, sample != 'anc1')

p1 <- ggplot(d_t4_590, aes(substrate, od_cor, col = interaction(evolved), text = sample)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotly::ggplotly(p1, tooltip = c('x', 'y', 'text'))
# 33 and 45 look like outliers

d_t4_590 <- filter(d_t4_590,! sample %in% c('33', '45'))

# add column for ranked mean OD_cor
d_t4_590 <- group_by(d_t4_590, substrate) %>%
  mutate(mean_od = mean(od_cor)) %>%
  ungroup() %>%
  mutate(., rank = dense_rank(desc(mean_od)),
         test = ifelse(plate <= 8, 'fuck', 'this'))

# plot performance across wells, ranked by best performance
plot1 <- ggplot(filter(d_t4_590, mean_od > 0.1)) +
  geom_line(aes(forcats::fct_reorder(substrate, rank), od_cor, group = sample, col = test), alpha = 0.25) +
  stat_summary(aes(rank, od_cor, col = test, group = evolved), fun.y = mean, geom = 'line') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = c(0.9, 0.8),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('optical density') +
  xlab('substrate rank') +
  ggtitle('Substrate rank across populations') +
  guides(col = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ evolved)

ggsave(file.path(path_fig, 'performance_plot.pdf'), plot1, height = 10, width = 10)
ggsave(file.path(path_fig, 'performance_plot2.pdf'), plot1a + facet_wrap(~evolved) + scale_color_manual(values = rep('black', 3)), height = 5, width = 12)


# Calculate phenotypic variance ####

# Take all the distance matrices of each pairwise combination of clones in each evolved
# average Euclidean distance
pop_dists_df <- filter(d_t4_590, evolved %in% c('with_community', 'without_community', 'ancestor')) %>%
  group_by(., evolved) %>%
  pairwise_dist(sample, substrate, od_cor, upper = FALSE) %>%
  rename(clone_i = item1, clone_j = item2)

# create average phenotypic diversity per population
# this is similar to a PCA and PCoA and betadisper()?
V_P <- group_by(pop_dists_df, evolved) %>%
  summarise(., V_P = mean(distance)) %>%
  data.frame()

# only two points here!!!

# calculate V_G - genotypic variance ####
# variance of OD across each genotype averaged over all the environments
V_G <- group_by(d_t4_590, evolved, substrate) %>%
  summarise(V_G = var(od_cor)) %>%
  data.frame()
V_G_pop <- group_by(V_G, evolved) %>%
  summarise(V_G = mean(V_G)) %>%
  data.frame()

# calculate V_E - environmental variance
# average variance of each clone across all the environments
V_E <- group_by(d_t4_590, evolved, sample) %>%
  summarise(V_E = var(od_cor)) %>%
  data.frame()
V_E_pop <- group_by(d_t4_590, evolved) %>%
  summarise(V_E = mean(V_E)) %>%
  data.frame()

# analyses
# Genotypic variance
mod_vg <- lmer(V_G ~ evolved + (1|substrate), V_G)
emmeans::emmeans(mod_vg, pairwise ~ evolved)
# NOPE

# environmental variance
mod_pg <- lm(V_E ~ evolved, V_E)
lsmeans::lsmeans(mod_pg, pairwise ~ evolved)

# plot genotypic and environmental variance across evolveds ####
# plot V_G and V_E ####
V_G_plot <- ggplot(V_G, aes(evolved, V_G)) +
  geom_boxplot(aes(fill = evolved, col = evolved), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(evolved, V_G, col = evolved), shape = 21, fill ='white', position = position_jitter(width = 0.1)) +
  ylab('genotypic variance') +
  xlab('evolved') +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle(expression(Genotypic~variance~(V[G]))) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)

V_E_plot <- ggplot(V_E, aes(evolved, V_E)) +
  geom_boxplot(aes(col = evolved, fill = evolved), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(evolved, V_E, col = evolved), shape = 21, fill ='white', position = position_jitter(width = 0.2)) +
  ylab('environmental variance') +
  xlab('evolved') +
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
  geom_density_ridges2(aes(x = OD_cor, y = factor(rank), fill = evolved, col = evolved), alpha = 0.5, rel_min_height = 0.01) +
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
d_sd <- group_by(d_t4_590, evolved, population, sample) %>%
  summarise(., sd_E = sd(od_cor)) %>%
  data.frame()

# create 2 copies of this for merging later
sd_j_clone <- dplyr::rename(d_sd, clone_j = sample, sd_j = sd_E) 
sd_i_clone <- rename(d_sd, clone_i = sample, sd_i = sd_E)

# create every pairwise combination of 1:n (clones/genotypes) for each population
d_R <- group_by(d_sd, evolved, population) %>%
  do(data.frame(expand.grid(clone_j = .$sample, clone_i = .$sample))) %>%
  ungroup() %>%
  filter(., clone_j > clone_i) %>%
  merge(., sd_j_clone, by = c('clone_j', 'evolved', 'population')) %>%
  merge(., sd_i_clone, by = c('clone_i', 'evolved', 'population'))

2# calculate R for each pairwise combination
d_R <- group_by(d_R, evolved, population) %>%
  mutate(., R_comb = (sd_j - sd_i)^2/(length(unique(clone_i))*(length(unique(clone_i))-1))) %>%
  ungroup()

# calculate responsiveness for each population
# sum of all the pairwise combinations
d_R_pop <- group_by(d_R, evolved, population) %>%
  summarise(., R_pop = sum(R_comb)) %>%
  data.frame()

# Plot responsiveness
r_plot <- ggplot(d_R_pop, aes(evolved, R_pop)) +
  geom_point(aes(evolved, R_pop, col = evolved), size = 3) +
  ylab('responsiveness') +
  xlab('evolved') +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('(a) Responsiveness') +
  scale_colour_viridis(discrete = TRUE)

# not significantly different
# summary(lm(R_pop ~ evolved, d_R_pop))

# calculate inconsistency ####

# prep data for calculating correlations
d_pearson <- group_by(d_t4_590, evolved, population) %>%
  do(data.frame(pairwise_cor(sample, substrate, od_cor, upper = FALSE))) %>%
  data.frame()

# merge dataframe to responsiveness dataframe
d_Inconsist <- merge(d_R, d_pearson, by = c('evolved', 'pop', 'clone_j', 'clone_i')) %>%
  mutate(., i = (sd_j*sd_i*(1-pear_cor))/(length(unique(clone_i))*(length(unique(clone_i))-1)))
d_I_pop <- group_by(d_Inconsist, evolved, pop) %>%
  summarise(., I_pop = sum(i),
            pear_pop = mean(pear_cor)) %>%
  data.frame()

# plot inconsistency
I_plot <- ggplot(d_I_pop, aes(evolved, I_pop)) +
  geom_point(aes(evolved, I_pop, col = evolved), size = 3) +
  ylab('Inconsistency') +
  xlab('evolved') +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('(b) Inconsistency') +
  scale_color_viridis(discrete = TRUE)

p_V_by_G <- r_plot + I_plot

ggsave(file.path(path_fig, 'responsiveness.pdf'), p_V_by_G, height = 5, width = 10)

summary(lm(I_pop ~ treat, d_I_pop))
 
# try a pca ####
d <- unite(d, 'id2', c(evolved, id), sep = '_', remove = FALSE)

# data set ready for PCA
d2 <- filter(d, id != 50 & id != 49)
d_PCA <- d2 %>%
  select(., starts_with('X'))
row.names(d_PCA) <- d2$id2

# create matrix
Euclid_mat <- dist(d_PCA)

# get variables for PCA
d_vars <- select(d2, id, id2, evolved) %>%
  mutate(., evolved = as.factor(evolved))
row.names(d_vars) <- d2$id2
PCA <- prcomp(d_PCA)
biplot(PCA)

# quick and dirty beta disper model
mod_adonis <- vegan::adonis(d_PCA ~ evolved, d_vars) # yes they have different centroids
mctoolsr::calc_pairwise_permanovas(Euclid_mat, d_vars, 'evolved')

Euclid_mat <- dist(d_PCA)
mod <- vegan::betadisper(Euclid_mat, d_vars$evolved)
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
  ggtitle('PCoA across evolveds') +
  guides(col = guide_legend(ncol = 8))

ggsave(file.path(path_fig, 'PCoA_across_evolveds.pdf'), last_plot(), height = 5, width = 7)

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
