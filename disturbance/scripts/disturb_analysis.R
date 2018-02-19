# example zero-inflated poisson model ####

# plots
library(rstan)
library(bayesplot)
library(brms)
library(ggplot2)
library(rstanarm)
library(tidybayes)

# set up model
model_bf <- bf(log.bet ~ Age*Sex + AgeSq*Sex + as.factor(TtI) + Community.Size + (1|Year) + (1 | ID), zi ~ Age*Sex + AgeSq*Sex + as.factor(TtI) + Community.Size + (1|Year) + (1 | ID))

# run model
model_brms <- brm(model_bf, 
                   data = data2, 
                   family = zero_inflated_poisson,
                   chains = 1,
                   iter = 100)

# summary of model
summary(model_brms)