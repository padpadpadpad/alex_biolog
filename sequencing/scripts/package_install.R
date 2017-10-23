# install all the packages necessary to run everything
# this is most likely overkill as I took the list of packages that encompasses some scripts and function we are less likely to use, but lets install them all anyway

# install devtools if necessary
install.packages('devtools')

# install MicrobioUoE
devtools::install_github('padpadpadpad/MicrobioUoE')

# list packages ####
# cran
cran_packages  <-  c("knitr", 
                     "phyloseqGraphTest", 
                     "phyloseq", 
                     "shiny",
                     "miniUI", 
                     "caret", 
                     "pls", 
                     "e1071", 
                     "ggplot2", 
                     "randomForest",
                     "vegan", 
                     "plyr", 
                     "dplyr",
                     "nlme",
                     "reshape2",
                     'tidyr',
                     'magrittr',
                     'vegan',
                     "PMA",
                     "structSSI",
                     "ade4",
                     "igraph", 
                     "ggnetwork", 
                     "intergraph", 
                     "scales", 
                     "phangorn")

# github
github_packages <- c("jfukuyama/phyloseqGraphTest", 
                     'benjjneb/dada2', 
                     'slowkow/ggrepel')

# bioconductor
bioc_packages <- c("phyloseq", 
                   "genefilter", 
                   "impute", 
                   "DECIPHER")

# install all packages
MicrobioUoE::package_install_all(cran_packages = cran_packages, github_packages = github_packages, bioc_packages = bioc_packages)

# will give a list of the failed packages
# The only one for me was phyloseqGraphTest and PMA
devtools::install_github('jfukuyama/phyloseqGraphTest')
install.packages('PMA')

# bingo - re-run to check
MicrobioUoE::package_install_all(cran_packages = cran_packages, github_packages = github_packages, bioc_packages = bioc_packages)
# Huzzah

# change this test
