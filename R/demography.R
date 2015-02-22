reproduce <- function(data, counter, params){
  
  # Subset so we're only dealing with alive adults
  adult_indices <- which(data$type == "adult" & data$alive == TRUE)
    # Skip reproduction if there are no adults!
  if(length(adult_indices) == 0 ) return(data) 
  
  n_adults_current <- length(adult_indices)
  
    # Calc number of total seedlings and the row indices we'll be using
  n_seedlings <- n_adults_current * params$crop_size 
  seedling_indices <- (counter$plant + 1):(counter$plant + n_seedlings)
  
  if(length(seedling_indices) != n_seedlings){
    stop("Problem with n_seedlings and indices")}
  
    # Update data frame with information for the seedlings  
  data$id[seedling_indices]      <-  seedling_indices
  data$type[seedling_indices]    <-  "seedling"
  data$alive[seedling_indices]   <-  TRUE
  data$age[seedling_indices]     <-  0
  
    # Seed dispersal process - adds coordinates and mother_id   
    data <- disperseSeed(data, params, adult_indices, seedling_indices)
    data <- dispersePollen(data, params, adult_indices, seedling_indices)

  return(data)
}


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



# Determine which plants survive
survival <- function(data, params){
    # If it's the first generation, skip
  if(counter$step == 0) next
  
    # If out of bounds and boundary setting is 'unsuitable', plants die
      # Could maybe speed up by only looking at plants that are currently alive
  if(params$boundary_setting == "unsuitable"){
    data$alive[data$pos_x > params$x_max | data$pos_x < 0] <- FALSE
    data$alive[data$pos_y > params$y_max | data$pos_y < 0] <- FALSE 
  }
    
    
    # Make a vector of whether adults that are currently alive will die
  adult_fate  <- sample(c(TRUE, FALSE), 
                         sum(with(data, alive[type == "adult" & alive == TRUE]),
                             na.rm = TRUE), replace = TRUE, 
                        prob = c(params$adult_survival, 1 - params$adult_survival))
  
  data$alive[with(data, which(type == "adult" & alive == TRUE))] <- adult_fate
  
    ## Seedling survival
  seedling_fate  <- sample(c(TRUE, FALSE), 
                      sum(with(data, alive[type == "seedling" & alive == TRUE]),
                             na.rm = TRUE), replace = TRUE, 
                      prob = c(params$seedling_survival, 1 - params$seedling_survival))
  
  data$alive[with(data, which(type == "seedling" & alive == TRUE))] <- seedling_fate
 
  return(data)  
}

# Increase age by one step
increaseAge <- function(data){
  alive <- which(data$alive == TRUE)
  data$age[alive] <- data$age[alive] + 1
  return(data)
}

# Transition to higher age class
transitionType <- function(data, params){
  data$type[data$type == "seedling" & data$age >= params$age_at_adult] <- "adult"
  return(data)
}







