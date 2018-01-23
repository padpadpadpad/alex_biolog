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
  mutate(., n_clones = case_when(treatment %in% c('4_related_clones', '4_unrelated_clones') ~ 4,
                                 treatment %in% c('individual_clone', 'lacz_ancestor', 'wt_ancestor') ~ 1,
                                 treatment %in% c('evolved_with_community', 'evolved_without_community') ~ 24,
                                 TRUE ~ NA_real_)) %>%
  arrange(., SampleID) %>% 
  select(., -c(X, num_clones))

write.csv(meta, 'sequencing/data/metadata.csv', row.names = FALSE)
