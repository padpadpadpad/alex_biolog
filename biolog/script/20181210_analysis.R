#-----------------------------------#
# analysis of new biolog plate data #
#-----------------------------------#

# clear workspace
rm(list = ls())

# load packages ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggridges)
library(lme4)
library(patchwork)
library(widyr)
library(corrr)
library(MicrobioUoE)
library(dglm)
library(hglm)

# figure path
path_fig <- 'plots'

# source extra functions ####
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

# look at variation across replicates of the same clone
d_sameclone <- filter(d, sample == 'anc1')
ggplot(d_sameclone, aes(substrate, od_cor, col = interaction(set, plate), shape = as.factor(od_wave))) +
  geom_point() +
  facet_wrap(~tp + od_wave, labeller = labeller(.multi_line = FALSE)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# look at variation across different ancestral clones
d_ancest <- filter(d, evolved == 'ancestor')
ggplot(d_ancest, aes(substrate, od_cor, col = sample, shape = as.factor(od_wave))) +
  geom_point(position = position_dodge(width = 1)) +
  facet_wrap(~ tp + od_wave, labeller = labeller(.multi_line = FALSE)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# lets use T4 590
d_t4_590 <- filter(d, od_wave == 590 & tp == 'T4')

# remove anc1 because it is super strange - this is a bit naughty
# it does something very different to everything else
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
plot1 <- ggplot(filter(d_t4_590, mean_od > 0.05)) +
  geom_line(aes(forcats::fct_reorder(substrate, rank), od_cor, group = sample, col = evolved), alpha = 1) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  ylab(expression(OD[600])) +
  xlab('substrate') +
  scale_color_manual(values = c('orange', 'dark grey', 'black')) +
  ggtitle('(a) Resource-use of pre-adapted clones')

ggsave(file.path(path_fig, 'performance_plot.pdf'), plot1, height = 10, width = 10)

# filter so that only wells where growth could occur were considered
d_t4_590 <- filter(d_t4_590, mean_od > 0.05)
# gives us 18 substrates

#----------------------------------#
# Calculate phenotypic variance ####
#----------------------------------#

# Take all the distance matrices of each pairwise combination of clones in each evolved
# average Euclidean distance
pop_dists_df <- filter(d_t4_590, evolved %in% c('with_community', 'without_community', 'ancestor')) %>%
  group_by(., evolved, population) %>%
  pairwise_dist(sample, substrate, od_cor, upper = FALSE) %>%
  rename(clone_i = item1, clone_j = item2)

# create average phenotypic diversity per population
# this is similar to a PCA and PCoA and betadisper()?
V_P <- group_by(pop_dists_df, evolved, population) %>%
  summarise(., V_P = mean(distance)) %>%
  data.frame()

# only two points here!!!

#---------------------------------------#
# calculate V_G - genotypic variance ####
#---------------------------------------#

# variance of OD across each genotype averaged over all the environments
# again does this need to be done across replicates rather than clones?
V_G <- group_by(d_t4_590, evolved, substrate, population) %>%
  summarise(V_G = var(od_cor)) %>%
  data.frame()
V_G_pop <- group_by(V_G, evolved, population) %>%
  summarise(V_G = mean(V_G)) %>%
  data.frame()

# calculate V_E - environmental variance
# average variance of each clone across all the environments
# does this need to be done on a replicate basis?
V_E <- group_by(d_t4_590, evolved, sample, population) %>%
  summarise(V_E = var(od_cor)) %>%
  data.frame()
V_E_pop <- group_by(V_E, evolved, population) %>%
  summarise(V_E = mean(V_E)) %>%
  data.frame()

# analyses
# Genotypic variance
mod_vg <- lmer(V_G ~ evolved + (1|substrate), V_G)
emmeans::emmeans(mod_vg, pairwise ~ evolved)
# NOPE

# environmental variance
mod_pg <- lm(V_E ~ evolved, V_E)
emmeans::emmeans(mod_vg, pairwise ~ evolved)

# plot genotypic and environmental variance across treatments ####
# plot V_G and V_E ####
V_G_plot <- ggplot(V_G_pop, aes(evolved, V_G)) +
  geom_pretty_boxplot(aes(evolved, V_G), filter(V_G_pop, evolved != 'ancestor'), col = 'black', fill = 'black') +
  geom_point(aes(evolved, V_G), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(V_G_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_G), shape = 21, fill = 'white', size = 5, filter(V_G_pop, evolved == 'ancestor')) +
  ylab('genotypic variance') +
  xlab('evolved') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none') +
  ggtitle(expression(Genotypic~variance~(V[G])))

V_E_plot <- ggplot(V_E_pop, aes(evolved, V_E)) +
  geom_pretty_boxplot(aes(evolved, V_E, col = evolved, fill = evolved), filter(V_E_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_E, col = evolved), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(V_E_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_E), fill = 'white', size = 7, col = 'orange', filter(V_E_pop, evolved == 'ancestor')) +
  ylab('environmental variance') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 12, color = 'black')) +
  ggtitle(expression((b)~Environmental~variance~(V[E]))) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))


V_P_plot <- ggplot(V_P, aes(evolved, V_P)) +
  geom_pretty_boxplot(aes(evolved, V_P), filter(V_P, evolved != 'ancestor'), col = 'black', fill = 'black') +
  geom_point(aes(evolved, V_P), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(V_P, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_P), shape = 21, fill = 'white', size = 5, filter(V_P, evolved == 'ancestor')) +
  ylab('phenotypic variance') +
  xlab('evolved') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none') +
  ggtitle(expression(Phenotypic~variance~(V[P])))


# plot
plot2 <- V_P_plot + V_G_plot + V_E_plot

#ggsave(file.path(path_fig, 'VE_VP_VG.pdf'), plot2, height = 5, width = 18)

#----------------------------------------------------#
# Calculate G x E interaction for each population ####
#----------------------------------------------------#

# see Barrett et al. 2005 Am Nat and Venail et al. 2008 Nature
# 1. calculate responsiveness - indicates differences in the environmental variances and thus measures diversity of resource exploitation strategies (specialists and generalists)
# sum (sd_j - sd_i)^2/(2*n_genotypes(n_genotypes - 1))

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

# calculate R for each pairwise combination
d_R <- group_by(d_R, evolved, population) %>%
  mutate(., R_comb = (sd_j - sd_i)^2/(2*n())*(n()-1)) %>%
  ungroup()

# calculate responsiveness for each population
# sum of all the pairwise combinations
d_R_pop <- group_by(d_R, evolved, population) %>%
  summarise(., R_pop = sum(R_comb)) %>%
  data.frame()

# Plot responsiveness
r_plot <- ggplot(d_R_pop, aes(evolved, R_pop)) +
  geom_pretty_boxplot(aes(evolved, R_pop, col = evolved, fill = evolved), filter(d_R_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, R_pop, col = evolved), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(d_R_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, R_pop), col = 'orange', size = 7, filter(d_R_pop, evolved == 'ancestor')) +
  ylab('responsiveness') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 12, color = 'black')) +
  ggtitle(expression((c)~Responsiveness)) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))


# not significantly different
# summary(lm(R_pop ~ evolved, d_R_pop))

#----------------------------#
# calculate inconsistency ####
#----------------------------#

# prep data for calculating correlations
d_pearson <- group_by(d_t4_590, evolved, population) %>%
  nest() %>%
  mutate(., cor_col = purrr::map(data, pairwise_cor, item = sample, feature = substrate, value = od_cor, upper = FALSE)) %>%
  unnest(cor_col) %>%
  rename(clone_i = item1, clone_j = item2) 

# merge dataframe to responsiveness dataframe
d_inconsist <- merge(d_pearson, sd_i_clone, by = c('evolved', 'population', 'clone_i'), all.x = TRUE) %>%
  merge(., sd_j_clone, by = c('evolved', 'population', 'clone_j'), all.x = TRUE) %>%
  group_by(., evolved, population) %>%
  mutate(., i = (sd_j*sd_i*(1-correlation))/(n()*(n()-1))) %>%
  summarise(., I_pop = sum(i),
            pear_pop = mean(correlation)) %>%
  data.frame()

# plot inconsistency
I_plot <- ggplot(d_inconsist, aes(evolved, I_pop)) +
  geom_pretty_boxplot(aes(evolved, I_pop, col = evolved, fill = evolved), filter(d_inconsist, evolved != 'ancestor')) +
  geom_point(aes(evolved, I_pop, col = evolved), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(d_inconsist, evolved != 'ancestor')) +
  geom_point(aes(evolved, I_pop), col = 'orange', size = 7, filter(d_inconsist, evolved == 'ancestor')) +
  ylab('Inconsistency') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 12, color = 'black')) +
  ggtitle(expression((d)~Inconsistency)) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))


phenotype_plot <- plot1 + {V_E_plot + r_plot + I_plot} + plot_layout(nrow = 2, heights = c(0.6, 0.4))

# save plot, other ways are available
ggsave(file.path(path_fig, 'biolog.png'), phenotype_plot, height = 12, width = 14)
ggsave(file.path(path_fig, 'biolog.pdf'), phenotype_plot, height = 12, width = 14)


p_V_by_G <- r_plot + I_plot

# ggsave(file.path(path_fig, 'V_GE_interaction.pdf'), p_V_by_G, height = 6, width = 15)

fit_inconsist <- lm(log(I_pop) ~ evolved, filter(d_inconsist, evolved != 'ancestor'))

# plot of fitness

# load data
d <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names()

d <- filter(d, treatment %in% c('individual_clone'))

# plot
ggplot(d, aes(evolution, fitness)) +
  MicrobioUoE::geom_pretty_boxplot(aes(col = evolution, fill = evolution)) +
  geom_point(aes(col = evolution), fill = 'white', shape = 21, position  = position_jitter(width = 0.1), size = 4) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(size = 16, color = 'black'),
        legend.position = 'none') +
  xlab('') +
  ylab('Relative fitness') +
  scale_x_discrete(labels = c('pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black')) +
  facet_wrap(~preadapt_pop)
  
# save plot, other ways are available
ggsave(file.path(path_fig, 'rel_fitness.png'), last_plot(), height = 5, width = 6)
ggsave(file.path(path_fig, 'rel_fitness.pdf'), last_plot(), height = 5, width = 6)

# try a normal linear mixed model
model <- lme4::lmer(fitness ~ evolution + (1|preadapt_pop), d)
cAIC4::cAIC(model)
model2 <- lme4::lmer(fitness ~ 1 + (1|preadapt_pop), d)
anova(model, model2)

# try a double generalised linear model
model_dglm <- dglm(fitness ~ evolution, ~1, data = d)
model_dglm2 <- dglm(fitness ~ evolution + preadapt_pop, ~evolution, data = d)
summary(model_dglm2)

# try a random double generalised linear model
model1 <- hglm(fixed = fitness ~ evolution,
              random = ~1|preadapt_pop,
              family = gaussian(link = "identity"),
              disp = ~ evolution,
              data = d,
              calc.like = TRUE)
model2 <- hglm(fixed = fitness ~ evolution,
               random = ~1|preadapt_pop,
               family = gaussian(link = "identity"),
               disp = ~ 1,
               data = d,
               calc.like = TRUE)

summary(model1)
summary(model2)
# essentially I do not know enough to do this. Maybe a simpler way. Would be cool if I can show this though

d_ranef <- data.frame(ranef =  model1$ranef, 
                      treat = rep(c('without_comm', 'with_comm'), each = 6)) 

boxplot(ranef~treat, d_ranef)

d_ave <- group_by(d, preadapt_pop, evolution) %>%
  summarise(mean = mean(fitness),
            sd = sd(fitness)) %>%
  ungroup()

gather(d_ave, 'metric', 'value', c(mean, sd)) %>%
  ggplot(., aes(evolution, value)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white') +
  facet_wrap(~metric, scales = 'free_y')

mod_mean <- lm(mean ~ evolution, d_ave)
summary(mod_mean)
mod_sd <- lm(sd ~ evolution, d_ave)
summary(mod_sd)

# correct each point by subtracting by the mean of each one
d <- group_by(d, preadapt_pop) %>%
  mutate(., fitness_cor = fitness - mean(fitness)) %>%
  ungroup()

ggplot(d, aes(evolution, fitness_cor)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, position  = position_jitter(width = 0.1), size = 4) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(size = 16)) +
  xlab('') +
  ylab('Relative fitness') +
  scale_x_discrete(labels = c('pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc'))

bartlett.test(fitness_cor ~ evolution, data = d)
car::leveneTest(fitness_cor ~ evolution, data = d)
