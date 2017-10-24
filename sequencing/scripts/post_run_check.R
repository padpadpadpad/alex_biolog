# first follow up script ####

library(ggplot2)
library(phyloseq)

# remove objects from file
rm(list - ls())

# set working directory to output folder of the run
setwd('sequencing/data/output/20171023_15:21')

# load in tracking of reads
tracking <- readRDS('20171023_15:21_track_reads_through_stages.rds')
samps <- row.names(tracking)
tracking <- data.frame(tracking) %>%
  mutate(samps = samps) %>%
  gather(., 'stage', 'reads', c(input, filtered, denoised_and_merged, tabled, nonchim))

head(tracking)

ggplot(tracking, aes(forcats::fct_relevel(stage, c('input', 'filtered', 'denoised_and_merged', 'tabled', 'nonchim')), reads)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2, height = 0), shape = 21, fill = 'white') +
  ylab('Number of reads') +
  xlab('Sequencing stage')

# average number of reads per sample per stage
track_ave <- group_by(tracking, stage) %>%
  summarise(., mean = mean(reads),
            min = min(reads),
            max = max(reads),
            sd = sd(reads)) %>%
  ungroup() %>%
  arrange(., desc(mean))

# average proportion of reads kept
min(track_ave$mean)/max(track_ave$mean)

# This might seem a little low - I would need to read around if this is bad or not

# lets plot a simple graph

# load in phyloseq object
ps <- readRDS('20171023_15:21_ps.rds')

# this take fucking ages
plot_bar(ps_prop, fill = "Genus")
