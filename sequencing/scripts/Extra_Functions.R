# Extra functions used in the analysis

# plot quality profiles
plot_qual <- function(filename, data_Fwd, data_Rev, subsample = NULL, ...){
  if(is.null(subsample)){subsample = c(1:length(data_Fwd))}
  pdf(filename, ...)
  for(i in subsample){print(plotQualityProfile(data_Fwd[i]) + ggtitle("Fwd"))}
  list_Rev <- list()
  for(i in subsample){print(plotQualityProfile(data_Rev[i]) + ggtitle("Rev"))}
  dev.off()
}

# initial setup
raw_read_setup <- function(packages, raw_path, filt_path, plot_path = NULL, output_path = NULL, progress_path = NULL, ref_fasta = NULL, ref_fasta_spp = NULL, meta_data = NULL, fwd_error = NULL, rev_error = NULL, run_filter){
  
  # assign some things to global environment
  assign('raw_path', raw_path, envir = globalenv())
  assign('filt_path', filt_path, envir = globalenv())
  assign('run_filter', run_filter, envir = globalenv())
  
  # create folders if they are not present
  if(!is.null(plot_path)){suppressWarnings(dir.create(plot_path))}
  if(!is.null(filt_path)){suppressWarnings(dir.create(filt_path))}
  if(!is.null(progress_path)){suppressWarnings(dir.create(progress_path))}
  
  # assign extra bits and pieces
  if(!is.null(fwd_error)){assign('fwd_error', readRDS(fwd_error), envir = globalenv())}
  if(!is.null(rev_error)){assign('rev_error', readRDS(rev_error), envir = globalenv())}
  if(!is.null(output_path)) assign('output_path', output_path, envir = globalenv())
  if(!is.null(ref_fasta)) assign('ref_fasta', ref_fasta, envir = globalenv())
  if(!is.null(ref_fasta_spp)) assign('ref_fasta_spp', ref_fasta_spp, envir = globalenv())
  
  if(!is.null(progress_path)) {assign('time', format(Sys.time(), '%Y%m%d_%H:%M_'), envir = globalenv())
  file.create(paste(progress_path, '/', time, 'progress.txt', sep = ''))
  assign('progress_file', paste(progress_path, '/', format(Sys.time(), '%Y%m%d_%H:%M_'), 'progress.txt', sep = ''), envir = globalenv())
  writeLines(paste('Run started at ', Sys.time()), progress_file)}
  if(!is.null(plot_path)){
    if(!file_test("-d", file.path(plot_path, substr(time, 1, nchar(time) - 1)))) dir.create(file.path(plot_path, substr(time, 1, nchar(time) - 1)))
        assign('plot_path', file.path(plot_path, substr(time, 1, nchar(time) - 1)), env = globalenv())
  }
  lapply(packages, library, character.only = TRUE)
  if(!is.null(meta_data)){
    assign('meta', read.csv(meta_data, stringsAsFactors = FALSE), envir = globalenv())
  }
  
}