# analysis script

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# figure path
path_fig <- 'figs'

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

# plot performance across wells, ranked by best performance
plot1a <- group_by(d_stack, id) %>%
arrange(., desc(OD)) %>%
  mutate(., rank = 1:96) %>%
  ggplot(.) +
  geom_line(aes(rank, OD, group = id), alpha = 0.25) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = 'none') +
  ylab('optical density') +
  xlab('substrate rank') +
  ggtitle('Substrate rank across populations')

# plot performance across well, without ranking by best performance
plot1b <- ggplot(d_stack) +
  geom_line(aes(C_source, OD, group = id), alpha = 0.25) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = 'none') +
  ylab('optical density') +
  xlab('substrate') +
  ggtitle('Performance across substrates')

plot1 <- gridExtra::grid.arrange(plot1a, plot1b, ncol = 1)

ggsave(file.path(path_fig, 'performance_plot.pdf'), plot1, height = 10, width = 10)

