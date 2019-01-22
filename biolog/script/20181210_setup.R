# setup biolog plates

# load in packages ####
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(janitor)
library(readr)

# function for reading in plates
# can be done for meta data or raw data
read_plate <- function(x, type = c('raw_data', 'meta_data')){
  if(missing(type)){type <- 'raw_data'}
  temp <- readxl::read_excel(x) %>%
    janitor::clean_names() %>%
    dplyr::select(1:13) %>%
    tidyr::gather(., 'well', 'od', 2:13) %>%
    dplyr::mutate(., well = readr::parse_number(well),
                  well = paste(x_1, well, sep = '_'),
                  file = basename(tools::file_path_sans_ext(x))) %>%
    dplyr::select(file, well, od)
  if(type == 'raw_data'){temp <- dplyr::mutate(temp, time = file.mtime(x))}
  
  if(type == 'meta_data'){temp <- rename(temp, treatment = od)}
  return(temp)
}

stock_sol_vol2 <- function (stock_sol_conc, new_sol_conc, diluent_vol)
{
  return((new_sol_conc * diluent_vol)/(stock_sol_conc - new_sol_conc))
}

round_to <- function(x, y) {
  if((y - x %% y) <= x %% y) { x + (y - x %% y)}
  else { x - (x %% y)}
}

# look inoculum
d <- read_plate('biolog/data/20181203_biolog_od_prenorm.xlsx') %>%
  mutate(., cont = 0.035,
         od_cor = od - cont,
         col = parse_number(well),
         row = gsub('_.*', '', well)) %>%
  filter(od > 0.04)

# get stock sol vols
d <- mutate(d, vol_sol = stock_sol_vol2(od_cor, 0.0015, 6000)) %>%
  mutate(vol_sol2 = round_to(vol_sol, 10))

select(d, row, col, vol_sol2) %>%
  spread(col, vol_sol2) %>%
  write.csv('biolog/data/20181203_biolog_od_sol_vols.csv')

select(d, row, col, od_cor) %>%
  spread(col, od_cor) %>%
  write.csv('data/20181022_TPC/od_cor.csv')
