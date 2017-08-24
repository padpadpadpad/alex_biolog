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
          'gridExtra')

# check if they are already installed
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]

# install packages
if(length(new_pkgs)) install.packages(new_pkgs)

# GitHub packages ####
# github pack age to install
github_packages <- 'padpadpadpad/MicrobioUoE'

# check if they are already installed
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]

# install packages
if(length(new_pkgs)) devtools::install_github(new_pkgs)
