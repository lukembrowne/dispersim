
disperseSeed <- function(sim, adult_indices, seedling_indices){
  
  # Calc dispersal distances
  disp_dist <- rweibull(length(seedling_indices), 
                        scale = sim$params$seed_kernel_scale, 
                        shape = sim$params$seed_kernel_shape)
  
  # Calculate direction that seed is moving in radians
  direction <- runif(length(seedling_indices), min = 0, max = 2*pi)
  
  # Coordinates of seed - not sure if this calc is right - needs testing!!
  x_coord <- sin(direction) * disp_dist
  y_coord <- cos(direction) * disp_dist
  
  # Save coords of adults, repeated X times depending on crop size
  x_coord_ad <- rep(sim$data$pos_x[adult_indices], each = sim$params$crop_size)
  y_coord_ad <- rep(sim$data$pos_y[adult_indices], each = sim$params$crop_size)
  
  # Put coordinates of seedlings into world dataframe
  sim$data$pos_x[seedling_indices] <- x_coord + x_coord_ad
  sim$data$pos_y[seedling_indices] <- y_coord + y_coord_ad
  
  # Add in mother id and color
  sim$data$id_mother[seedling_indices] <- rep(sim$data$id[adult_indices], 
                                          each = sim$params$crop_size)
  sim$data$color[seedling_indices] <- rep(sim$data$color[adult_indices], 
                                      each = sim$params$crop_size)
}


dispersePollen <- function(sim, adult_indices, seedling_indices){
  
  # Calc distance for each mother to potential pollen donors
  dist_mat <-  as.matrix(dist(cbind(sim$data$pos_x[adult_indices], 
                                    sim$data$pos_y[adult_indices]), upper = T))
  diag(dist_mat) <- NA
  rownames(dist_mat) <- adult_indices
  colnames(dist_mat) <- adult_indices
  
  # Calculate probabilty of fertilization based on dispersal kernel
  # Choose pollen donor based on the probability of fertilization  
  prob_mat <- 1 - pweibull(dist_mat[], shape = sim$params$pollen_kernel_shape,
                       scale = sim$params$pollen_kernel_scale
                       )
  
  # Arranged as matrix with columns as mothers, rows as pollen donors
  pollen_donors <- apply(prob_mat, 1, function(x){
    x  <- x[!is.na(x)]
    # Sample for pollen donors based on probability
    as.numeric(sample(names(x), sim$params$crop_size, 
                      replace = TRUE, prob = x))
  })
  # Assign id_father to that seedling
  sim$data$id_father[seedling_indices] <- as.numeric(pollen_donors) # Squash matrix
}


## Calculate seed dispersal distance
# Distance from seedling to mother

calcDispDistance <- function(sim, type, method, step){
  
  if(type != "adult" & type != "seedling"){
    stop("Type not recognized...\n")
  }
  
  if(method != "seed" & method != "pollen" & method != "total pollen"){
    stop("Method not recognized...\n")
  }
  
    # Subset to relevant plants... those that are alive at time step of interest
  registry_index <- ifelse(type == "adult", 1, 2)
  alive <- sim$data[sim$registry[[step]][[registry_index]], ]
  
    # Distance from plant to mother
  if(method == "seed"){
    x1 <- alive$pos_x
    y1 <- alive$pos_y
    x2 <- sim$data$pos_x[alive$id_mother]
    y2 <- sim$data$pos_y[alive$id_mother]
  }
  
    # Distance from mother to father
  if(method == "pollen"){
    x1 <- sim$data$pos_x[alive$id_father]
    y1 <- sim$data$pos_y[alive$id_father]
    x2 <- sim$data$pos_x[alive$id_mother]
    y2 <- sim$data$pos_y[alive$id_mother]
  }
  
    # Total pollen movement - from plant to father
  if(method == "total pollen"){
    x1 <- alive$pos_x
    y1 <- alive$pos_y
    x2 <- sim$data$pos_x[alive$id_father]
    y2 <- sim$data$pos_y[alive$id_father]
  }
  
  disp_dist <- .calcDist(x1 = x1, 
                       x2 = x2,
                       y1 = y1, 
                       y2 = y2)
  
  return(disp_dist)  # Returns vector of distances 
}
