# raw read processing ####

# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# clean workspace before starting ####
rm(list = ls())

# set seed ####
set.seed(42)

# start time
start_time <- Sys.time()

# setup paths and packages and error if previously estimated ####
# this will automatically set up the environment to run the analysis
MicrobioUoE::dada2_raw_read_setup(meta_data = 'data/metadata_example.csv')

# get time
time <- paste0(basename(plot_path), '_', collapse = '')

cat(paste('\nThis run is done using raw_read_processing.R'), file = progress_file, append = TRUE)

# list files ####
fns <- sort(list.files(raw_path, pattern = 'fast', full.names = TRUE, recursive = T))

# sort files for forward and reverse sequences ####
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# name the files in the list ####
# need to be the same as in the metadata file
sample_namesF <- paste('sample', gsub('_.*', '', basename(fnFs)), sep = '_')
sample_namesR <- paste('sample', gsub('_.*', '', basename(fnRs)), sep = '_')

# check Fwd and Rev files are in the same order
if(!identical(sample_namesF, sample_namesR)) stop("Forward and reverse files do not match.")

# check Fwd and Rev files are in the SampleID column of the meta data
if(all(sample_namesF %in% meta$SampleID) == FALSE) stop("Forward and reverse file names are not present in metadata column SampleID")

# run filter parameters ####
# this can be based on the quality profiles in qual_plot_preFilt.pdf
if(run_filter == 'Y'){
  
  # create filenames
  mast_pre_filtF <- paste0(file.path(output_path, 'fwd_master_preFilt_'), basename(plot_path), '.fastq', collapse = '')
  mast_pre_filtR <- paste0(file.path(output_path, 'rev_master_preFilt_'), basename(plot_path), '.fastq', collapse = '')
  
  # create master forward and reverse files ####
  system(paste('cat', paste(fnFs, collapse = ' '), '>', mast_pre_filtF, sep = ' '))
  system(paste('cat', paste(fnRs, collapse = ' '), '>', mast_pre_filtR, sep = ' '))
  
  # check quality of data ####
  pdf(file.path(plot_path, 'qual_plot_preFilt.pdf'))
  print(plotQualityProfile(mast_pre_filtF, n = 2e6) +
    ggtitle('Fwd reads master quality profile'))
  print(plotQualityProfile(mast_pre_filtR, n = 2e6) +
    ggtitle('Rev reads master quality profile'))
  dev.off()
  
  # Trim and filter ####
  filtFs <- file.path(filt_path, basename(fnFs)) 
  filtRs <- file.path(filt_path, basename(fnRs))
  
  # set trimming parameters ####
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       trimLeft = 25, 
                       truncLen=c(250, 250),
                       maxN = 0, 
                       maxEE = 2,
                       truncQ = 2,
                       compress=TRUE,
                       multithread = TRUE)
  
  # plot quality post filtering ####
  mast_post_filtF <- paste0(file.path(output_path, 'fwd_master_postFilt_'), basename(plot_path), '.fastq', collapse = '')
  mast_post_filtR <- paste0(file.path(output_path, 'rev_master_postFilt_'), basename(plot_path), '.fastq', collapse = '')
  
  # create master forward and reverse files
  system(paste('cat', paste(filtFs, collapse = ' '), '>', mast_post_filtF, sep = ' '))
  system(paste('cat', paste(filtRs, collapse = ' '), '>', mast_post_filtR, sep = ' '))
  
  pdf(file.path(plot_path, 'qual_plot_postFilt.pdf'))
  print(plotQualityProfile(mast_post_filtF, n = 2e6) +
    ggtitle('Fwd reads master quality profile'))
  print(plotQualityProfile(mast_post_filtR, n = 2e6) +
    ggtitle('Rev reads master quality profile'))
  dev.off()
  
  # add update to progress file
  cat(paste('\nFiltering completed at ', format(Sys.time(), '%Y-%m-%d %H:%m')), file = progress_file, append = TRUE)
  
  # remove some objects
  rm(list = c('mast_pre_filtF', 'mast_post_filtF', 'mast_pre_filtR', 'mast_post_filtR'))
  
}

# check the number of reads post filtering
head(out)

# get filt path ####
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# name files
names(filtFs) <- sample_namesF
names(filtRs) <- sample_namesR

# from here this should be run in one go
# learn error rates ####
# if this has been done before, assign error rates
if(length(grep('fwd_error', ls())) == 1){
  errF <- fwd_error
  errR <- rev_error
}

# if this has not been done before, infer error rates
if(length(grep('fwd_error', ls())) == 0){
  
  # learn error rates ####
  # learn forward error rates
  dd_learnF <- learnErrors(filtFs, 
                           multithread = TRUE,
                           randomize = TRUE,
                           MAX_CONSIST = 30)
  
  cat(paste('\nForward error rates completed at', Sys.time()), file = progress_file, append = TRUE)
  cat(paste('\nForward error rate:', dada2:::checkConvergence(dd_learnF), sep = ' '), file = progress_file, append = TRUE)
  # Learn reverse error rates
  dd_learnR <- learnErrors(filtRs, 
                           multithread = TRUE,
                           randomize = TRUE,
                           MAX_CONSIST = 30)
  
  cat(paste('\nReverse error rates completed at', Sys.time()), file = progress_file, append = TRUE)
  cat(paste('\nReverse error rates:', dada2:::checkConvergence(dd_learnR), sep = ' '), file = progress_file, append = TRUE)
  
  # save out error rates
  saveRDS(dd_learnF, paste(output_path, '/', time, 'fwd_error.rds', sep = ''))
  saveRDS(dd_learnR, paste(output_path, '/', time, 'rev_error.rds', sep = ''))
}

# dereplicate files
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_namesF
names(derepRs) <- sample_namesR
cat(paste('\nDereplicated forward and reverse files', Sys.time()), file = progress_file, append = TRUE)

# infer sequence data ####
dadaFs <- dada(derepFs, err = dd_learnF, pool = TRUE, multithread = TRUE)
cat(paste('\nForward sequences inferred', Sys.time()), file = progress_file, append = TRUE)
dadaRs <- dada(derepRs, err = dd_learnR, pool = TRUE, multithread = TRUE)
cat(paste('\nReverse sequences inferred', Sys.time()), file = progress_file, append = TRUE)

# save out the dada class objects
saveRDS(dadaFs, paste(output_path, '/', time, 'dadaFs.rds', sep = ''))
saveRDS(dadaRs, paste(output_path, '/', time, 'dadaRs.rds', sep = ''))

# merge inferred forward and reverse sequences ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
cat(paste('\nInferred forward and reverse sequences merged', Sys.time()), file = progress_file, append = TRUE)

# construct sequence table ####
seqtab_all <- makeSequenceTable(mergers)
cat(paste('\nSequence table constructed', Sys.time()), file = progress_file, append = TRUE)

# remove chimeric sequences ####
seqtab <- removeBimeraDenovo(seqtab_all)
cat(paste('\nChimeric sequences removed', Sys.time()), file = progress_file, append = TRUE)

# track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab_all), rowSums(seqtab))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample_namesF
head(track)

saveRDS(track, paste(output_path, '/', time, 'track_reads_through_stages.rds', sep = ''))

# plot out tracking of sample reads through stages ####
samps <- row.names(track)
track <- data.frame(track) %>%
  mutate(samps = samps) %>%
  gather(., 'stage', 'reads', c(input, filtered, denoised, merged, tabled, nonchim))

ggplot(track, aes(forcats::fct_relevel(stage, c('input', 'filtered', 'denoised', 'merged', 'tabled', 'nonchim')), reads)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2, height = 0), shape = 21, fill = 'white') +
  ylab('Number of reads') +
  xlab('Sequencing stage')

ggsave(paste(plot_path, '/', time, 'track_reads.pdf', sep = ''), height = 5, width = 7)


# assign taxonomy ####
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
cat(paste('\nTaxonomy assigned', Sys.time()), file = progress_file, append = TRUE)
if(!is.null(ref_fasta_spp)){
  cat(paste('\nAssigning species at', Sys.time(), '\n'), file = progress_file, append = TRUE)
  spp_assign <- capture.output(taxtab <- addSpecies(taxtab, refFasta = ref_fasta_spp, verbose = TRUE))
  cat(spp_assign, file = progress_file, append = TRUE)
  cat(paste('\nSpecies assigned', Sys.time()), file = progress_file, append = TRUE)
}
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# save files in case phylogeny does not run
saveRDS(taxtab, paste(output_path, '/', time, 'taxtab.rds', sep = ''))
saveRDS(seqtab, paste(output_path, '/', time, 'seqtab.rds', sep = ''))

# subset meta for just the samples present

# SampleID needs to match sample_namesF
meta <- filter(meta, SampleID %in% sample_namesF)
rownames(meta) <- meta$SampleID
ps <- phyloseq(tax_table(taxtab), 
               sample_data(meta),
               otu_table(seqtab, taxa_are_rows = FALSE))

# save a phyloseq object without the phylogeny
saveRDS(ps, paste(output_path, '/', time, 'ps.rds', sep = ''))
save(ps, file = paste(output_path, '/', time, 'ps.Rdata', sep = ''))

# plot error rates ####
pdf(paste(plot_path, '/', time, 'error_rates.pdf', sep = ''))
plotErrors(dd_learnF) +
  ggtitle('Forward error rates')
plotErrors(dd_learnR) +
  ggtitle('Reverse error rates')
dev.off()

cat(paste('\nEnd of raw read processing without construction of phylogeny', Sys.time()), file = progress_file, append = TRUE)

# End time
end_time <- Sys.time()

cat(paste('\nThis run (without phylogeny estimation) took:', difftime(end_time, start_time, unit = 'hours'), 'hours', sep = ' '), file = progress_file, append = TRUE)

# multiple alignment ####
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
cat(paste('\nSequences aligned', Sys.time()), file = progress_file, append = TRUE)

# construct phylogenetic tree using phangorn ####
phang_align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA") 
dm <- phangorn::dist.ml(phang_align)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
fit <-  phangorn::pml(treeNJ, data = phang_align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                              rearrangement = "stochastic", control = phangorn::pml.control(trace = 0))
cat(paste('\nConstructed phylogenetic tree', Sys.time()), file = progress_file, append = TRUE)

# save files
saveRDS(taxtab, paste(output_path, '/', time, 'taxtab.rds', sep = ''))
saveRDS(seqtab, paste(output_path, '/', time, 'seqtab.rds', sep = ''))
saveRDS(phang_align, paste(output_path, '/', time, 'phytree.rds', sep = ''))
# files saved

# subset meta for just the samples present
meta <- filter(meta, SampleID %in% sample_namesF)
rownames(meta) <- meta$SampleID
ps <- phyloseq(tax_table(taxtab), 
               sample_data(meta),
               otu_table(seqtab, taxa_are_rows = FALSE), 
               phy_tree(fitGTR$tree))

# save phyloseq object
saveRDS(ps, paste(output_path, '/', time, 'ps.rds', sep = ''))
save(ps, file = paste(output_path, '/', time, 'ps.Rdata', sep = ''))

cat(paste('\nEnd of raw read processing', Sys.time()), file = progress_file, append = TRUE)

# End time
end_time <- Sys.time()

cat(paste('\nThis run took:', difftime(end_time, start_time, unit = 'hours'), 'hours', sep = ' '), file = progress_file, append = TRUE)