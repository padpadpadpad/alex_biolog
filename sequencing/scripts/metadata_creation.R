# create metadata file ####

# list files

# currently is only a column with sample_1 etc in for Sample_ID - this needs to be expanded

meta <- data.frame(SampleID = paste('sample', 1:102, sep = '_'), stringsAsFactors = FALSE)

write.csv(meta, 'sequencing/data/metadata.csv', row.names = FALSE)