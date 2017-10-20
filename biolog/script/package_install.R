# install all necessary packages 

# first need to install devtools to install the MicrobioUoE package

# CRAN packages ####
# package list in scripts/analysis.R
pkgs <- c('tidyr',
          'dplyr',
          'ggplot2',
          'nlme',
          'readxl',
          'devtools',
          'gridExtra',
          'viridis',
          'ggridges')

# check if they are already installed
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]

# install packages
if(length(new_pkgs)) install.packages(new_pkgs)

# GitHub packages ####
# github packages to install
github_packages <- 'padpadpadpad/MicrobioUoE'

# check if they are already installed
new_pkgs <- github_packages[!(github_packages %in% installed.packages()[,"Package"])]

# install packages
if(length(new_pkgs)) devtools::install_github(new_pkgs)
