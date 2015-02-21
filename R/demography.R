reproduce <- function(data, counter, crop_size){
  
  # Subset so we're only dealing with alive adults
  adult_indices <- which(data$type == "adult" & data$alive == TRUE)
    # Skip reproduction if there are no adults!
  if(length(adult_indices) == 0 ) return(data) 
  
  n_adults_current <- length(adult_indices)
  
    # Calc number of total seedlings and the row indices we'll be using
  n_seedlings <- n_adults_current * crop_size 
  seedling_indices <- (counter$plant + 1):(counter$plant + n_seedlings)
  
  if(length(seedling_indices) != n_seedlings){
    stop("Problem with n_seedlings and indices")}
  
    # Update data frame with information for the seedlings  
  data$id[seedling_indices]      <-  seedling_indices
  data$type[seedling_indices]    <-  "seedling"
  data$alive[seedling_indices]   <-  TRUE
  data$age[seedling_indices]     <-  0
  
    # Seed dispersal process - adds coordinates and mother_id   
    data <- disperseSeed(data, adult_indices, seedling_indices, crop_size)
  
  # Increase plant counter
  # Could add test that last row without NA should also equal where the plant counter is
  counter$plant <<- counter$plant + n_seedlings
  
  return(data)
}


disperseSeed <- function(data, adult_indices, seedling_indices, crop_size){
  
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
  x_coord_ad <- rep(data$pos_x[adult_indices], each = crop_size)
  y_coord_ad <- rep(data$pos_y[adult_indices], each = crop_size)

    # Put coordinates of seed into world dataframe
  data$pos_x[seedling_indices] <- x_coord + x_coord_ad
  data$pos_y[seedling_indices] <- y_coord + y_coord_ad
  
    # Add in mother id
  data$id_mother[seedling_indices] <- rep(data$id[adult_indices], each = crop_size)
  data$color[seedling_indices] <- rep(data$color[adult_indices], each = crop_size)
  
  return(data)   
  
}


survival <- function(data){
    # If it's the first generation, skip
  if(counter$step == 0) next
    
    # Make a vector of whether adults that are currently alive will die
  adult_fate  <- sample(c(TRUE, FALSE), 
                         sum(with(data, alive[type == "adult" & alive == TRUE]),
                             na.rm = TRUE), replace = TRUE, prob = c(.9, .1))
  data$alive[with(data, which(type == "adult" & alive == TRUE))] <- adult_fate
  
  seedling_fate  <- sample(c(TRUE, FALSE), 
                      sum(with(data, alive[type == "seedling" & alive == TRUE]),
                             na.rm = TRUE), replace = TRUE, prob = c(.5, .5))
  data$alive[with(data, which(type == "seedling" & alive == TRUE))] <- seedling_fate
 
  return(data)  
}




#survival