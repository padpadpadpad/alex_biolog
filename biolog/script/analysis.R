# analysis script

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggridges)

# figure path
path_fig <- 'biolog/figs'

# source extra functions
source('biolog/script/functions.R')

# load in data ####
d <- MicrobioUoE::bind_biolog_all('biolog/data/biolog_data_assay_1.xlsx', sheets = 'Sheet1')

# make into long format ####

# meta data
meta <- data.frame(id = 1:50, treatment = c(rep('comm', times = 24), rep('no_comm', times = 24), rep('wild_type', times = 2)), stringsAsFactors = FALSE)

# data and metadata together
d <- merge(d, meta, by = 'id')

# load in ancestral data
d_ancest <- MicrobioUoE::bind_biolog_all('biolog/data/20170124_Ancestors_gn2biolog.xlsx', sheets = 'Sheet1') %>%
  mutate(id = id + 50)
meta_ancest <- data.frame(id = 51:58, treatment = 'wild_type', stringsAsFactors = FALSE)
d_ancest <- merge(d_ancest, meta_ancest, by = 'id')
d <- bind_rows(d, d_ancest)

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
  stat_summary(aes(rank, OD_cor, col = treatment), fun.y = mean, geom = 'line', lwd = 1.25) +
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
  ylim(c(0, 0.1))

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

# plot all carbon sources
ggplot(d_stack2) +
  geom_density_ridges2(aes(x = OD_cor, y = factor(rank), fill = treatment, col = treatment), alpha = 0.5, rel_min_height = 0.01) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_point(aes(x = WT, y = factor(rank)), size = 0.5) +
  theme_bw()

ggsave(file.path(path_fig, 'crazy_ggjoy_plot.pdf'), last_plot(), height = 12, width = 6)

# try a pca
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
  theme(legend.position = c(.83, .85)) +
  xlab('') +
  scale_x_discrete(labels = c('Community', 'No Community', 'Wild Type'))

mod <- lm(distances ~ group, betadisper_dat$distances)
summary(mod)   
