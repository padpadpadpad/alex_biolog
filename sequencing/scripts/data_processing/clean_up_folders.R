# clean up folders and progress files ####

# when we are happy with our runs, we have also likely created many folders in the data/output, data/progress and plots folders that could be reduced

# in this script I list all the folders in these systematically, and delete them if they are under a certain size, indicating that the full workflow did not work

# load packages ####
library(dplyr)
library(tidyr)
library(ggplot2)

# files to keep here
files_to_keep <- c('data/progress/progress_folder.txt',
                   'data/output/output_folder.txt',
                   'plots/plots_folder.txt',
                   'data/progress/example_progress_big_data.txt',
                   'data/progress/example_progress_raw_reads.txt')

# progress files ####
progress_files <- data.frame(file = list.files('data/progress', full.names = TRUE), stringsAsFactors = FALSE) %>%
  mutate(., size = file.size(file)) %>%
  filter(! file %in% files_to_keep) %>%
  filter(size < 1000)

# look at progress files. Generally delete files below 1 kb
folder_size <- function(files){
  temp <- sum(file.size(list.files(files, all.files = TRUE, full.names = TRUE, recursive = TRUE)))
  return(temp)
}

# output folder ####
output_folder <- data.frame(file = list.files('data/output', full.names = TRUE), stringsAsFactors = FALSE) %>%
  filter(! file %in% files_to_keep) %>%
  rowwise() %>%
  mutate(., size = folder_size(file))

# look at output folder - set memory under which to delete folder
output_to_delete <- 3e5

output_folder <- filter(output_folder, size < output_to_delete)

# plots folder ####
plots_folder <- data.frame(file = list.files('plots', full.names = TRUE), stringsAsFactors = FALSE) %>%
  filter(! file %in% files_to_keep) %>%
  rowwise() %>%
  mutate(., size = folder_size(file))

# look at plots folder - set memory under which to delete folder
plots_to_delete <- 100000

plots_folder <- filter(plots_folder, size < plots_to_delete)

# combine all files and folders to delete
to_delete <- c(plots_folder$file, output_folder$file, progress_files$file)

# delete said files - WARNING THIS CANNOT BE UNDONE. RUN WITH CAUTION
unlink(to_delete, recursive = TRUE, force = TRUE)
