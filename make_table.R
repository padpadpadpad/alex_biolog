# create table 1

library(flextable)
library(htmltools)
library(readxl)

table <- read_excel('~/Google Drive/work/pseudomonas_diversity_preadaptation/figures/table_1.xlsx') %>%
  janitor::clean_names() %>%
  flextable() %>%
  autofit() %>%
  merge_v(., 'diversity_level') %>%
  theme_box()

htmltools_value(table)

html_print(table)
