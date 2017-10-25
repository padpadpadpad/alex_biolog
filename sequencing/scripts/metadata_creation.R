# create metadata file ####

# SampleID needs to match names of samples in big_data_processing or raw_read_processing
meta <- read.csv('sequencing/data/soil_community_data_samplenames.csv', stringsAsFactors = FALSE) %>%
  dplyr::rename(SampleID = Sample_number, 
         sample_name = Sample_name, 
         treatment = Treatment_label) %>%
  dplyr::mutate(treatment = tolower(treatment),
                treatment = gsub('-', '', treatment),
                SampleID = paste('sample', SampleID, sep = '_')) %>%
  arrange(., SampleID)

write.csv(meta, 'sequencing/data/metadata.csv', row.names = FALSE)