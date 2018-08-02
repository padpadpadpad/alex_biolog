# permute actual distance between samples from output from

# get_dist
get_dist <- function(d_group, d_centroid){
  dist <- dist(rbind(dplyr::filter(d_centroid, group == d_group$var1) %>% dplyr::pull(value), dplyr::filter(d_centroid, group == d_group$var2) %>% dplyr::pull(value)))
  return(dist)
}

all_combs <- function(group){
  d_group <- t(combn(unique(group), 2)) %>%
    data.frame(stringsAsFactors = FALSE)
  colnames(d_group) <- c('var1', 'var2')
  return(d_group)
}

# permute distances
permute_distance <- function(dist_mat, group, nperm = 9999){
  nperm = nperm
  
  temp <- tibble::tibble(dist = rep(NA, nperm + 1), type = rep(NA, nperm + 1))
  
  # do betadisper and extract centroids
  mod <- betadisper(dist_mat, group)
  d_centroid <- betadisper_centroids(mod) %>%
    gather(., 'axis', 'value', starts_with('PCoA'))
  
  # get dataframe with all combinations of group
  d_group <- all_combs(d_centroid$group)
  
  # calculate distance 
  act_dist <- get_dist(d_group, d_centroid)
  
  temp[1,1] <- act_dist
  temp[1,2] <- 'actual_distance'
  
  # bootstrap samples
  # positions of each sample
  pos <- which(group %in% d_group$var1)
  pos2 <- 1:length(group)
  pos2 <- pos2[!pos2 %in% pos]
  
  # calculate indexing of distance matrix
  ind_mat<-matrix(d_group$var1, nr = nrow(as.matrix(dist_mat)),nc = ncol(as.matrix(dist_mat)))
  for(i in 1:nrow(ind_mat)){
    for(j in 1:ncol(ind_mat)){
      if(i %in% pos2 == TRUE & j %in% pos2 == TRUE){ind_mat[i,j] <- d_group$var2}
    }
  }
  ind_mat[lower.tri(ind_mat)]
  
  for(i in 1:nperm){
    # create distance matrix
    dist_mat2 <- as.matrix(dist_mat)
    dist_mat2[upper.tri(dist_mat2, diag = TRUE)] <- NA
    dist_mat2[pos, pos] <- sample(dist_mat2[colnames(dist_mat2) %in% vec1 & rownames(dist_mat2) %in% vec1], replace = TRUE)
    dist_mat2[colnames(dist_mat2) %in% vec2 & rownames(dist_mat2) %in% vec2] <- sample(dist_mat2[colnames(dist_mat2) %in% vec2 & rownames(dist_mat2) %in% vec2], replace = TRUE)
    dist_mat <- as.dist(dist_mat2)

    
    
    temp[i+1, 1] <- act_dist2
  }
  
  temp[2:(nperm+1), 2] <- 'permuted'
  
  return(temp)
  
}

ind.mat<-matrix(0,nr=nrow(as.matrix(dist_mat)),nc=ncol(as.matrix(dist_mat)))

for(i in 1:nrow(ind.mat)){
  for(j in 1:ncol(ind.mat)){
    if(i%in%pos==TRUE&j%in%pos==TRUE){ind.mat[i,j]<-1}
  }
}

length(ind.mat[lower.tri(ind.mat)])

x <- ind.mat[lower.tri(ind.mat)]
