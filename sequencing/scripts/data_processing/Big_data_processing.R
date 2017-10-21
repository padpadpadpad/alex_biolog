# raw read Big Data ####

# setup ####
rm(list = ls())

# set seed ####
set.seed(42)

# start time
start_time <- Sys.time()

# load source script with additional functions
source('DNA/scripts/Extra_Functions.R')

# setup paths and packages
raw_read_setup(
  packages = c('ggplot2', 'gridExtra', 'dada2', 'phyloseq', 'DECIPHER', 'tidyr', 'dplyr'),
  raw_path = 'DNA/data/Trimmed',
  filt_path = 'DNA/data/Filtered',
  plot_path = 'DNA/plots',
  output_path = 'DNA/data/output',
  progress_path = 'DNA/data/progress',
  ref_fasta = "DNA/data/train_sets/clustering/final_combine/ref_db_final.fasta",
  meta_data = "DNA/data/pond_transplant_DNA_metadata.csv",
  fwd_error = "DNA/data/error/20161210_14:26_errorF.rds",
  rev_error = "DNA/data/error/20161210_14:26_errorR.rds",
  run_filter = 'N'
)

fns <- sort(list.files(raw_path, pattern = 'fast', full.names = TRUE, recursive = T))
# sort files for forward and reverse sequences ####
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# run filter parameters ####
if(run_filter == 'Y'){
  
  # check quality of data ####
  # subsample <- sample(length(fnFs), 3)
  if(!file.exists(file.path(plot_path, 'qual_plot_preFilt.pdf'))){
    plot_qual(file.path(plot_path, 'qual_plot_preFilt.pdf'), fnFs, fnRs, height = 5, width = 7)
  }
  
  # Trim and filter ####
  # create filter path if it is not present
  if(!file_test("-d", filt_path)) dir.create(filt_path) 
  filtFs <- file.path(filt_path, basename(fnFs)) 
  filtRs <- file.path(filt_path, basename(fnRs))
  
  # set trimming parameters ####
  for(i in seq_along(fnFs)) {
    fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]), 
                      c(filtFs[[i]], filtRs[[i]]),
                      trimLeft=0, truncLen=c(250, 250),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE)
  }
  
  plot_qual(file.path(plot_path, 'qual_plot_postFilt.pdf'), filtFs, filtRs, height = 5, width = 7)
  
  cat(paste('\nFiltering completed at ', format(Sys.time(), '%Y-%m-%d %H:%m')), file = progress_file, append = TRUE)
  
}

# get filt path
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# name the objects in the list
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(filtFs) <- sample.names
names(filtRs) <- sample.names

if(length(grep('fwd_error', ls())) == 1){
  errF <- fwd_error
  errR <- rev_error
}

if(length(grep('fwd_error', ls())) == 0){
  
  # learn error rates ####
  # learn forward error rates
  NSAM.LEARN <- 10 # Choose enough samples to have at least 1M reads
  drp.learnF <- derepFastq(sample(filtFs, NSAM.LEARN))
  dd.learnF <- dada(drp.learnF, err=NULL, selfConsist=TRUE, multithread=TRUE)
  errF <- dd.learnF[[1]]$err_out
  cat(paste('\nForward error rates completed at', Sys.time()), file = progress_file, append = TRUE)
  cat(paste('\nForward error rate:', dada2:::checkConvergence(dd.learnF[[1]]), sep = ' '), file = progress_file, append = TRUE)
  # Learn reverse error rates
  drp.learnR <- derepFastq(sample(filtRs, NSAM.LEARN))
  dd.learnR <- dada(drp.learnR, err=NULL, selfConsist=TRUE, multithread=TRUE)
  errR <- dd.learnR[[1]]$err_out
  cat(paste('\nReverse error rates completed at', Sys.time()), file = progress_file, append = TRUE)
  cat(paste('\nReverse error rates:', dada2:::checkConvergence(dd.learnR[[1]]), sep = ' '), file = progress_file, append = TRUE)

  # plot error rates ####
  pdf(paste(plot_path, '/', time, 'error_rates.pdf', sep = ''))
  plotErrors(dd.learnF)
  plotErrors(dd.learnR)
  dev.off()

  rm(drp.learnF);rm(dd.learnF)
  rm(drp.learnR);rm(dd.learnR)
}


# Sample inference and merger of paired-end reads ####
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(i in 1:length(sample.names)) {
  cat(paste("\nProcessing:",  i, 'of', length(sample.names), ':', sample.names[i], Sys.time(), sep = ' '), file = progress_file, append = TRUE)
  derepF <- derepFastq(filtFs[[sample.names[i]]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sample.names[i]]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sample.names[i]]] <- merger
}

rm(derepF); rm(derepR)

# construct sequence table ####
seqtab.all <- makeSequenceTable(mergers)
cat(paste('\nSequence table constructed', Sys.time()), file = progress_file, append = TRUE)

# remove chimeric sequences ####
seqtab <- removeBimeraDenovo(seqtab.all)
cat(paste('\nChimeric sequences removed', Sys.time()), file = progress_file, append = TRUE)

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
phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA") 
dm <- phangorn::dist.ml(phang.align)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
fit <-  phangorn::pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                              rearrangement = "stochastic", control = phangorn::pml.control(trace = 0))
cat(paste('\nConstructed phylogenetic tree', Sys.time()), file = progress_file, append = TRUE)

# save files
if(!file_test("-d", file.path(output_path, substr(time, 1, nchar(time) - 1)))) dir.create(file.path(output_path, substr(time, 1, nchar(time) - 1)))
output_path <- file.path(output_path, substr(time, 1, nchar(time) - 1))
saveRDS(taxtab, paste(output_path, '/', time, 'taxtab.rds', sep = ''))
saveRDS(taxtab, paste(output_path, '/', time, 'seqtab.rds', sep = ''))
saveRDS(taxtab, paste(output_path, '/', time, 'phytree.rds', sep = ''))
# files saved
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


