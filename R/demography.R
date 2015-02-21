reproduce <- function(world, crop_size = world$crop_size){
  
  # Subset so we're only dealing with alive adults
  adult_indices <- which(world$data$type == "adult" & world$data$alive == TRUE)
  if(length(adult_indices) == 0 ) return(world) # Skip reproduction if there are no adults!
  
  n_adults_current <- length(adult_indices)
  
    # Calc number of total seedlings and the row indices we'll be using
  n_seedlings <- n_adults_current * crop_size 
  seedling_indices <- (world$plant_counter + 1):(world$plant_counter + n_seedlings)
  
  if(length(seedling_indices) != n_seedlings) stop("Problem with n_seedlings and indices")
  
    # Update world data frame with information for the seedlings  
  world$data$id[seedling_indices]      <-  seedling_indices
  world$data$type[seedling_indices]    <-  "seedling"
  world$data$alive[seedling_indices]   <-  TRUE
  world$data$age[seedling_indices]     <-  0
  
    # Seed dispersal process - adds coordinates and mother_id
  world <- disperseSeed(world, adult_indices, seedling_indices, crop_size)
  
  # Increase plant counter
  # Could add test that last row without NA should also equal where the plant counter is
  world$plant_counter <- world$plant_counter + n_seedlings
  
  return(world)
}


disperseSeed <- function(world, adult_indices, seedling_indices, crop_size){
  
    # Calc dispersal distances
  disp_dist <- rweibull(length(seedling_indices), 
                        scale = world$seed_kernel_scale, 
                        shape = world$seed_kernel_shape)
  
    # Calculate direction that seed is moving in radians
  direction <- runif(length(seedling_indices), min = 0, max = 2*pi)
  
    # Coordinates of seed - not sure if this calc is right - needs testing!!
  x_coord <- sin(direction) * disp_dist
  y_coord <- cos(direction) * disp_dist
    
    # Save coords of adults 
  x_coord_ad <- rep(world$data$pos_x[adult_indices], each = crop_size)
  y_coord_ad <- rep(world$data$pos_y[adult_indices], each = crop_size)

    # Put coordinates of seed into world dataframe
  world$data$pos_x[seedling_indices] <- x_coord + x_coord_ad
  world$data$pos_y[seedling_indices] <- y_coord + y_coord_ad
  
    # Add in mother id
  world$data$id_mother[seedling_indices] <- rep(world$data$id[adult_indices],
                                                each = crop_size)
  world$data$color[seedling_indices] <- rep(world$data$color[adult_indices],
                                            each = crop_size)
  
  return(world)   
  
}


survival <- function(world){
    # If it's the first generation, skip
  if(world$step_counter == 0) next
    
    # Make a vector of whether adults that are currently alive will die
  adult_fate  <- sample(c(TRUE, FALSE), 
                         sum(with(world$data, alive[type == "adult" & alive == TRUE]),
                             na.rm = TRUE), replace = TRUE, prob = c(.9, .1))
  world$data$alive[with(world$data,
                        which(type == "adult" & alive == TRUE))] <- adult_fate
  
  seedling_fate  <- sample(c(TRUE, FALSE), 
                      sum(with(world$data, alive[type == "seedling" & alive == TRUE]),
                             na.rm = TRUE), replace = TRUE, prob = c(.5, .5))
  world$data$alive[with(world$data,
                        which(type == "seedling" & alive == TRUE))] <- seedling_fate
  
  
  
  return(world)
  
  
}




#survival