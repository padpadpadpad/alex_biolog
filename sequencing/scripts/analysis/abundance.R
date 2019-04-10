# look at metadata across treatments

# load packages
library(ggplot2)
library(dplyr)

# figure path
path_fig <- 'plots'

# load data
d <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names()

d_sub <- filter(d, !treatment %in% c('wt_ancestor', 'nmc_t0', 'negative_control')) %>%
  mutate(., log_density = log10(density_cfu_g)) %>%
  mutate_all(., function(x)ifelse(x == 0|is.infinite(x)|is.nan(x), NA, x))

# plot
ggplot(d_sub, aes(treatment, log_density)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, position  = position_jitter(width = 0.1), size = 3) +
  theme_bw() +
  ggtitle('Abundance of P. fluorescens at the end of the experiment') +
  xlab('Treatment') +
  ylab('Density (CFU / g)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(path_fig, 'pseudomonas_abundance.png'), last_plot(), height = 6, width = 8)

mod1 <- lm(log10(density_cfu_g) ~ treatment, d_sub)

# pairwise of pre_adaptation treatments
d_preadapt <- filter(d_sub, !treatment %in% c('4_unrelated_clones', 'lacz_ancestor')) %>%
  select(log_density, evolution, n_clones) %>%
  filter(!is.na(log_density))

ggplot(d_preadapt, aes(evolution, log_density)) +
  geom_point() +
  facet_wrap(~ n_clones)

mod_preadapt <- d_preadapt %>%
  nest(-n_clones) %>%
  mutate(., model = purrr::map(data, ~lm(log_density ~ evolution, .x))) %>%
  unnest(model %>% purrr::map(broom::tidy)) %>%
  filter(term == 'evolutionwithout_community') %>%
  mutate(., p.value = round(p.value, 3))

mod_preadapt

# pool all these together to look at the other effect
d_nclones <- select(d_sub, log_density, treatment, n_clones) %>%
  mutate(., n_clones = ifelse(treatment == 'lacz_ancestor', 'lacz_ancestor', n_clones)) %>%
  filter(!is.na(log_density))

ggplot(d_nclones, aes(n_clones, log_density)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, size = 3)

# problem with doing unbalanced anova is usually homogeneity of variances assumption
group_by(d_nclones, n_clones) %>%
  summarise(var = var(log_density),
            sd = sd(log_density))
# even though lacZ has lower n, has higher variance than individual clone where sample size is much bigger
# so can continue with normal lm and normal anovas

mod_treat <- lm(log_density ~ n_clones, d_nclones)
mod_treat2 <- lm(log_density ~ 1, d_nclones)
anova(mod_treat, mod_treat2)
emmeans::emmeans(mod_treat, pairwise ~ n_clones)

# create plot
ggplot(d_nclones, aes(forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')), log_density, col = forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')), fill = forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')))) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab(expression(log[10]~Density~(cfu~g^-1~soil))) +
  xlab('') +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones', 'negative\ncontrol')) +
  ylim(c(4, 7.5))

ggsave(file.path(path_fig, 'abundance.png'), last_plot(), height = 5, width = 6)
ggsave(file.path(path_fig, 'abundance.pdf'), last_plot(), height = 5, width = 6)

