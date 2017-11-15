# create metadata file ####

# load in packages ####
library(dplyr)
library(tidyr)

# SampleID needs to match names of samples in big_data_processing or raw_read_processing
meta <- read.csv('sequencing/data/soil_community_data_samplenames.csv', stringsAsFactors = FALSE) %>%
  dplyr::rename(SampleID = Sample_number, 
         sample_name = Sample_name, 
         treatment = Treatment_label,
         evolution = Evolved_with.without_community,
         preadapt_pop = Preadaptation.population,
         fitness = Fitness.W.,
         bio_var = Biolog_Variance) %>%
  dplyr::mutate(treatment = tolower(treatment),
                treatment = gsub('-', '', treatment),
                SampleID = paste('sample', SampleID, sep = '_')) %>%
  dplyr::mutate_at(., c('evolution', 'preadapt_pop'), function(x){ifelse(x == '', NA, x)}) %>%
  arrange(., SampleID) %>% 
  select(., -X)

write.csv(meta, 'sequencing/data/metadata.csv', row.names = FALSE)
