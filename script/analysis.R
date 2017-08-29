# analysis script

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggjoy)

# figure path
path_fig <- 'figs'

# source extra functions
source('script/functions.R')

# load in data ####
d <- MicrobioUoE::bind_biolog_all('data/biolog_data_assay_1.xlsx', sheets = 'Sheet1')

# make into long format ####

# meta data
meta <- data.frame(id = 1:50, treatment = c(rep('comm', times = 24), rep('no_comm', times = 24), rep('wild_type', times = 2)), stringsAsFactors = FALSE)

# data and metadata together
d <- merge(d, meta, by = 'id')

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

# plot performance across wells, ranked by best performance
plot1a <- group_by(d_stack, id) %>%
arrange(., desc(OD_cor)) %>%
  mutate(., rank = 1:96) %>%
  ggplot(.) +
  geom_line(aes(rank, OD_cor, group = id, col = treatment), alpha = 0.75) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = c(0.9, 0.8)) +
  ylab('optical density') +
  xlab('substrate rank') +
  ggtitle('Substrate rank across populations') +
  guides(col = guide_legend(override.aes = list(alpha = 1)))

# plot performance across well, without ranking by best performance
plot1b <- ggplot(d_stack) +
  geom_line(aes(C_source, OD_cor, group = id, col = treatment), alpha = 0.4) +
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
V_G_plot <- ggplot(filter(V_G, treatment != 'wild_type'), aes(treatment, V_G)) +
  geom_boxplot(aes(fill = treatment, col = treatment), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(treatment, V_G, col = treatment), shape = 21, fill ='white', position = position_jitter(width = 0.1)) +
  ylab('genotypic variance') +
  xlab('Treatment') +
  scale_x_discrete(labels = c('Community', 'No Community')) +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle(expression(Genotypic~variance~(V[G]))) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)

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
gridExtra::grid.arrange(V_G_plot, V_E_plot, ncol = 2)

# plot all carbon sources
ggplot(d_noWT) +
  geom_joy2(aes(x = OD_cor, y = factor(rank), fill = treatment, col = treatment), alpha = 0.5, rel_min_height = 0.01) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_point(aes(x = WT, y = factor(rank)), size = 0.5) +
  theme_bw()
