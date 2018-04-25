# create multi page pdf documents

# load staplr
library(staplr)

staple_pdf <- function(input_directory = NULL, input_files = NULL, output_filename = "Full_pdf", 
                        output_directory = NULL) 
{
  if(is.null(input_directory) & is.null(input_files)) {
    input_directory <- tcltk::tk_choose.dir(caption = "Select directory which contains PDF fies")
  }
  if(!is.null(input_directory)){input_filepaths <- (Sys.glob(file.path(input_directory, "*.pdf")))}
  if(!is.null(input_files)){input_filepaths <- input_files}
  
  if (is.null(output_directory)) {
    output_directory <- tcltk::tk_choose.dir(caption = "Select directory to save output")
  }
  output_filepath <- file.path(output_directory, paste(output_filename, 
                                                       ".pdf", sep = ""))
  quoted_names <- paste0("\"", input_filepaths, "\"")
  file_list <- paste(quoted_names, collapse = " ")
  output_filepath <- paste0("\"", output_filepath, "\"")
  system_command <- paste("pdftk", file_list, "cat", "output", 
                          output_filepath, sep = " ")
  system(command = system_command)
}

# list files to staple
files <- c('biolog/figs/old_vs_new_biologs.pdf', 
           'biolog/figs/old_vs_new_biologs_control.pdf',
           'biolog/figs/geno_enviro_var_plot.pdf')

# bind into one pdf
staple_pdf(input_files = files, output_directory = 'biolog/figs')
