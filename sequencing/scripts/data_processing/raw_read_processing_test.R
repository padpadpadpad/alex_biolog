# raw processing ####

# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# clean workspace before starting ####
rm(list = ls())

# set seed ####
set.seed(42)

# start time
start_time <- Sys.time()

# load source script with additional functions ####
source('sequencing/scripts/Extra_Functions.R')

# setup paths and packages and error if previously estimated ####
# this will automatically create a progress file each time it is run in the progress path
# will create relevant folders if they are not there - raw_path needs to be present as it has the data
raw_read_setup(
  packages = c('ggplot2', 'gridExtra', 'dada2', 'phyloseq', 'DECIPHER', 'tidyr', 'dplyr'),
  raw_path = 'sequencing/data/raw_data',
  filt_path = 'sequencing/data/Filtered',
  plot_path = 'sequencing/plots',
  output_path = 'sequencing/data/output',
  progress_path = 'sequencing/data/progress',
  ref_fasta = 'sequencing/data/ref_trainsets/rdp_train_set_16.fa',
  meta_data = 'sequencing/data/metadata.csv',
  fwd_error = NULL,
  rev_error = NULL,
  run_filter = 'Y'
)

# create folder for output
dir.create(file.path(output_path, substr(time, 1, nchar(time) - 1)))
output_path <- file.path(output_path, substr(time, 1, nchar(time) - 1))

# list files ####
fns <- sort(list.files(raw_path, pattern = 'fast', full.names = TRUE, recursive = T))

# sort files for forward and reverse sequences ####
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# for test pick the first 15 of each of these
fnFs <- fnFs[1:5]
fnRs <- fnRs[1:5]

# check quality of data ####
if(!file.exists(file.path(plot_path, 'qual_plot_preFilt_test.pdf'))){
  plot_qual(file.path(plot_path, 'qual_plot_preFilt_test.pdf'), fnFs, fnRs, height = 5, width = 7)
}

# run filter parameters ####
# this can be based on the quality profiles in qual_plot_preFilt.pdf
if(run_filter == 'Y'){
  
  # Trim and filter ####
  filtFs <- file.path(filt_path, basename(fnFs)) 
  filtRs <- file.path(filt_path, basename(fnRs))
  
  # set trimming parameters ####
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       trimLeft=0, truncLen=c(250, 250),
                       maxN=0, maxEE=2, truncQ=2,
                       compress=TRUE,
                       multithread = TRUE)
  
  # plot quality post filtering
  plot_qual(file.path(plot_path, 'qual_plot_postFilt_test.pdf'), filtFs, filtRs, height = 5, width = 7)
  
  # add update to progress file
  cat(paste('\nFiltering completed at ', format(Sys.time(), '%Y-%m-%d %H:%m')), file = progress_file, append = TRUE)
  
}

# check how many reads made it through the filtering step
head(out)

# get filt path ####
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# name the objects in the list
sample.names <- paste('sample', gsub('_.*', '', basename(filtFs)), sep = '_')
names(filtFs) <- sample.names
names(filtRs) <- sample.names

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
                           MAX_CONSIST = 2)
  errF <- dd_learnF$err_out
  cat(paste('\nForward error rates completed at', Sys.time()), file = progress_file, append = TRUE)
  cat(paste('\nForward error rate:', dada2:::checkConvergence(dd_learnF), sep = ' '), file = progress_file, append = TRUE)
  # Learn reverse error rates
  dd_learnR <- learnErrors(filtRs, 
                           multithread = TRUE,
                           randomize = TRUE,
                           MAX_CONSIST = 2)
  errR <- dd_learnR$err_out
  cat(paste('\nReverse error rates completed at', Sys.time()), file = progress_file, append = TRUE)
  cat(paste('\nReverse error rates:', dada2:::checkConvergence(dd_learnR), sep = ' '), file = progress_file, append = TRUE)
  
  dev.off()
  # plot error rates ####
  pdf(paste(plot_path, '/', time, 'error_rates.pdf', sep = ''))
  plotErrors(dd_learnF) +
    ggtitle('Forward error rates')
  plotErrors(dd_learnR) +
    ggtitle('Reverse error rates')
  dev.off()
  
  # save out error rates
  saveRDS(dd_learnF, paste(output_path, '/', time, 'fwd_error.rds', sep = ''))
  saveRDS(dd_learnR, paste(output_path, '/', time, 'rev_error.rds', sep = ''))
}

# dereplicate files
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
cat(paste('\nDereplicated forward and reverse files', Sys.time()), file = progress_file, append = TRUE)

# infer sequence data ####
dadaFs <- dada(derepFs, err = errF, pool = TRUE, multithread = TRUE)
cat(paste('\nForward sequences inferred', Sys.time()), file = progress_file, append = TRUE)
dadaRs <- dada(derepRs, err = errR, pool = TRUE, multithread = TRUE)
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
rownames(track) <- sample.names
head(track)

saveRDS(track, paste(output_path, '/', time, 'track_reads_through_stages.rds', sep = ''))

# assign taxonomy ####
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
cat(paste('\nTaxonomy assigned', Sys.time()), file = progress_file, append = TRUE)

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
saveRDS(taxtab, paste(output_path, '/', time, 'seqtab.rds', sep = ''))
saveRDS(taxtab, paste(output_path, '/', time, 'phytree.rds', sep = ''))
# files saved

# subset meta for just the samples present
meta <- filter(meta, SampleID %in% sample.names)
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

rm(list = setdiff(ls(), c("ps", 'taxtab', 'meta', 'seqtab', 'fitGTR')))
