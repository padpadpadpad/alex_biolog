# look at metadata across treatments

# load packages
library(ggplot2)
library(dplyr)

# figure path
path_fig <- 'sequencing/plots/fresh'

# load data
d <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names()

d <- filter(d, !treatment %in% c('wt_ancestor', 'nmc_t0'))

# plot
ggplot(d, aes(treatment, density_cfu_g)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, position  = position_jitter(width = 0.1), size = 3) +
  theme_bw() +
  ggtitle('Abundance of P. fluorescens at the end of the experiment') +
  xlab('Treatment') +
  ylab('Density (CFU / g)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "l")

ggsave(file.path(path_fig, 'pseudomonas_abundance.png'), last_plot(), height = 6, width = 8)

# look at variation in fitness and abundance

d_single <- filter(d, treatment == 'individual_clone')

ggplot(d_single, aes(fitness, density_cfu_g, col = evolution)) +
  geom_point()
