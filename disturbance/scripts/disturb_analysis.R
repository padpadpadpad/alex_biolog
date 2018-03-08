# example zero-inflated poisson model ####

# plots
library(rstan)
#library(bayesplot)
library(brms)
#library(rstanarm)
#library(tidybayes) # devtools::install_github('mjskay/tidybayes')
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(vegan)

# load in data
d <- read.csv('disturbance/data/Disturbance_data.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  mutate(., disturb_past = readr::parse_number(pre_treatment),
         rep = sub("[^[:alpha:]]+", "", pre_treatment)) %>%
  rename(disturb_current = post_treatment) %>%
  select(., -c(pre_treatment, pre_treatment2)) %>%
  select(., rep, disturb_past, disturb_current, everything()) %>%
  mutate(., same_disturb = ifelse(disturb_past == disturb_current, 'yes', 'no'),
         pres_abs = ifelse(post_competitor > 0, 1, 0))
d <- gather(d, morph, count, c(starts_with('pre_'))) %>%
  group_by(rep, disturb_past) %>%
  summarise(., inv_simp = diversity(count, index = 'invsimpson'),
            simp = diversity(count, index = 'simpson'),
            shannon = diversity(count, index = 'shannon')) %>%
  ungroup() %>%
  merge(., d, by = c('rep', 'disturb_past'))

gather(d, metric, diversity, c(inv_simp, simp, shannon)) %>%
  ggplot(., aes(as.factor(disturb_past), diversity)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = 'free')

d_tots <- 
  group_by(d, disturb_past, disturb_current, rep, same_disturb) %>%
  summarise(tot_before = pre_ws + pre_sm + pre_fs,
            tot_end = post_ws + post_sm + post_fs)  %>%
  ungroup()

d_tots2 <- gather(d_tots, time, count, starts_with('tot_'))
d_comp <- select(d, disturb_current, disturb_past, post_competitor, rep) %>%
  mutate(time = 'tot_end')

ggplot(d_tots2, aes(time, count, col = same_disturb)) +
  geom_point() +
  geom_line(aes(group = interaction(disturb_current, disturb_past, rep))) +
  geom_point(aes(time, post_competitor), col = 'black', d_comp) +
  facet_wrap(~ disturb_current)
  
# plot 
ggplot(d, aes(disturb_past, log10(post_competitor), col = same_disturb)) +
  geom_point() +
  facet_wrap(~ disturb_current, ncol = 5)

# try normal glmm
fit_glm <- glm(post_competitor ~ disturb_past*disturb_current, family = 'quasipoisson', d)


# set up model
model_bf <- bf(post_competitor ~ disturb_past*disturb_current, zi ~ disturb_past*disturb_current)

# run model
model_brms <- brm(model_bf, 
                   data = d, 
                   family = zero_inflated_poisson,
                   chains = 3,
                   iter = 1000)

# summary of model
summary(model_brms)