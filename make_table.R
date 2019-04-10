# create table 1

# load packages
library(flextable)
library(readxl)
library(rmarkdown)
library(webshot)

# create table
table <- read_excel('~/Google Drive/work/pseudomonas_diversity_preadaptation/figures/table_1.xlsx') %>%
  mutate_all(as.character) %>%
  janitor::clean_names() %>%
  flextable() %>%
  merge_v(., 'diversity_level') %>%
  merge_at(., i = c(1,2), j = 'number_of_replicates') %>%
  merge_at(., i = c(3,4), j = 'number_of_replicates') %>%
  merge_at(., i = c(6,7), j = 'number_of_replicates') %>%
  merge_v(., 'treatment_3') %>%
  align(j = c('diversity_level', 'number_of_replicates'), align = 'center', part = 'all') %>%
  set_header_labels(diversity_level = 'Diversity level',
                    treatment_2 = 'Treatment',
                    number_of_replicates = 'Number of replicates') %>%
  merge_h_range(., i = c(1,2), j1 = 'treatment_2', j2 = 'treatment_3', part = 'body') %>%
  merge_h_range(., i = c(6:9), j1 = 'treatment_2', j2 = 'treatment_3', part = 'body') %>%
  merge_h_range(., j1 = 'treatment_2', j2 = 'treatment_3', part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  align(j = c('treatment_2', 'treatment_3'), align = 'left', part = 'body') %>%
  align(j = c('treatment_2'), align = 'center', part = 'header') %>%
  hline(i = c(2, 4, 5, 7, 8), border = officer::fp_border(color="black"))

# save as a png
# create an Rmd file
rmd_name <- tempfile(fileext = ".Rmd")
cat("```{r echo=FALSE}\ntable\n```", file = rmd_name)

# render as an html file ----
html_name <- tempfile(fileext = ".html")
render(rmd_name, output_format = "html_document", output_file = html_name )

# get a png from the html file with webshot ----
webshot(html_name, zoom = 2, file = "table_1.png", 
        selector = "body > div.container-fluid.main-container > div.tabwid > table")
