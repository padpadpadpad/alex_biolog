# extra functions to use in the analysis

# function for getting Euclidean distances of each population
pop_dists <- function(x){
  temp <- dist(dplyr::select(x, starts_with('X')))
  temp_mat <- reshape2::melt(as.matrix(temp), varnames = c('row', 'col')) %>% 
    dplyr::filter(., row > col)
  return(temp_mat)
}

# function for getting environmental correlation of each pair of genotypes
pop_cor <- function(x, n_gen, rows_delete){
  temp <- data.frame(t(x), stringsAsFactors = FALSE) %>%
    dplyr::filter(., ! row.names(.) %in% rows_delete) %>%
    dplyr::mutate_all(., as.numeric)
  colnames(temp) <- 1:n_gen
  # delete rows
  temp_mat <- as.matrix(cor(temp)) %>%
    reshape2::melt(as.matrix(.), varnames = c('clone_j', 'clone_i')) %>%
    dplyr::filter(., clone_j > clone_i) %>%
    rename(., pear_cor = value)
  return(temp_mat)
}

# get distance from 00
dist_from_00 <- function(x, y){
  return(sqrt((0 - x)^2+(0-y)^2))
}

# getting distances from betadisper() object
betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}

# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# betadisper data
get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}

