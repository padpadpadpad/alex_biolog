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
