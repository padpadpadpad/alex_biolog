# run error rate estimation and save values
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
  ref_fasta = "DNA/data/train_sets/phytoRef_trainset.fasta",
  meta_data = "DNA/data/pond_transplant_DNA_metadata.csv",
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

# learn error rates ####
# learn forward error rates
NSAM.LEARN <- length(filtFs) # Choose enough samples to have at least 1M reads
drp.learnF <- derepFastq(sample(filtFs, NSAM.LEARN))
dd.learnF <- dada(drp.learnF, err=NULL, selfConsist=TRUE, multithread=TRUE, MAX_CONSIST = 25)
errF <- dd.learnF[[1]]$err_out
cat(paste('\nForward error rates completed at', Sys.time()), file = progress_file, append = TRUE)
cat(paste('\nForward error rate:', dada2:::checkConvergence(dd.learnF[[1]]), sep = ' '), file = progress_file, append = TRUE)
# Learn reverse error rates
drp.learnR <- derepFastq(sample(filtRs, NSAM.LEARN))
dd.learnR <- dada(drp.learnR, err=NULL, selfConsist=TRUE, multithread=TRUE,MAX_CONSIST = 25)
errR <- dd.learnR[[1]]$err_out
cat(paste('\nReverse error rates completed at', Sys.time()), file = progress_file, append = TRUE)
cat(paste('\nReverse error rates:', dada2:::checkConvergence(dd.learnR[[1]]), sep = ' '), file = progress_file, append = TRUE)

# save output
saveRDS(errR, paste(output_path, '/', time, 'errorR.rds', sep = ''))
saveRDS(errF, paste(output_path, '/', time, 'errorF.rds', sep = ''))
