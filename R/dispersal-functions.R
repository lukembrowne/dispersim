
disperseSeed <- function(data, params, adult_indices, seedling_indices){
  
  # Calc dispersal distances
  disp_dist <- rweibull(length(seedling_indices), 
                        scale = params$seed_kernel_scale, 
                        shape = params$seed_kernel_shape)
  
  # Calculate direction that seed is moving in radians
  direction <- runif(length(seedling_indices), min = 0, max = 2*pi)
  
  # Coordinates of seed - not sure if this calc is right - needs testing!!
  x_coord <- sin(direction) * disp_dist
  y_coord <- cos(direction) * disp_dist
  
  # Save coords of adults 
  x_coord_ad <- rep(data$pos_x[adult_indices], each = params$crop_size)
  y_coord_ad <- rep(data$pos_y[adult_indices], each = params$crop_size)
  
  # Put coordinates of seed into world dataframe
  data$pos_x[seedling_indices] <- x_coord + x_coord_ad
  data$pos_y[seedling_indices] <- y_coord + y_coord_ad
  
  # Add in mother id
  data$id_mother[seedling_indices] <- rep(data$id[adult_indices], 
                                          each = params$crop_size)
  data$color[seedling_indices] <- rep(data$color[adult_indices], 
                                      each = params$crop_size)
  
  return(data)     
}


dispersePollen <- function(data, params, adult_indices, seedling_indices){
  
  # Calc distance for each mother to potential pollen donors
  dist_mat <-  as.matrix(dist(cbind(data$pos_x[adult_indices], 
                                    data$pos_y[adult_indices]), upper = T))
  diag(dist_mat) <- NA
  rownames(dist_mat) <- adult_indices
  colnames(dist_mat) <- adult_indices
  
  # Calculate probabilty of fertilization based on dispersal kernel
  # Choose pollen donor based on the probability of fertilization  
  prob_mat <- pweibull(dist_mat[], shape = params$pollen_kernel_shape,
                       scale = params$pollen_kernel_scale)
  # Arranged as matrix with columns as mothers, rows as pollen donors
  pollen_donors <- apply(prob_mat, 1, function(x){
    x  <- x[!is.na(x)]
    as.numeric(sample(names(x), params$crop_size, 
                      replace = TRUE, prob = x))
  })
  # Assign id_father to that seedling
  data$id_father[seedling_indices] <- as.numeric(pollen_donors) # Squash matrix
  return(data)
}









# Calculate distance between two points
calcDist <- function(x1, x2, y1, y2){
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}

## Calculate seed dispersal distance
# Distance from seedling to mother

calcSeedDispDistance <- function(data){
  
  recruits <- data[!is.na(data$id_mother), ]
  
  disp_dist <- calcDist(x1 = recruits$pos_x, x2 = data$pos_x[recruits$id_mother],
                        y1 = recruits$pos_y, y2 = data$pos_y[recruits$id_mother])
  
  return(disp_dist)  # Returns vector of distances 
}

# Calculate pollen dispersal distance
# Distance from mother to father
calcPollenDispDistance <- function(data){
  
  recruits <- data[!is.na(data$id_father), ]
  
  disp_dist <- calcDist(x1 = data$pos_x[recruits$id_father], 
                        x2 = data$pos_x[recruits$id_mother],
                        y1 = data$pos_y[recruits$id_father], 
                        y2 = data$pos_y[recruits$id_mother])
  
  return(disp_dist)  # Returns vector of distances 
}



# Calculate 'effective' pollen dispersal distance -
# Distance from seedling to father
calcEffectivePollenDispDistance <- function(data){
  
  recruits <- data[!is.na(data$id_father), ]
  
  disp_dist <- calcDist(x1 = recruits$pos_x, 
                        x2 = data$pos_x[recruits$id_father],
                        y1 = recruits$pos_y, 
                        y2 = data$pos_y[recruits$id_father])
  
  return(disp_dist)  # Returns vector of distances 
}

