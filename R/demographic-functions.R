procReproduce <- function(sim){
  
  # Subset so we're only dealing with alive adults
  adult_indices <- which(sim$data$type == "adult" & sim$data$alive == TRUE)
    # Skip reproduction if there are no adults!
  if(length(adult_indices) == 0 ) return(cat("No adults alive left to breed...\n")) 
  
  n_adults_current <- length(adult_indices)
  
    # Calc number of total seedlings and the row indices we'll be using
  n_seedlings <- n_adults_current * sim$params$crop_size 
  seedling_indices <- (sim$counter$plant + 1):(sim$counter$plant + n_seedlings)
  
  if(length(seedling_indices) != n_seedlings){
    stop("Problem with n_seedlings and indices")}
  
    # Update data frame with information for the seedlings  
  sim$data$id[seedling_indices]      <-  seedling_indices
  sim$data$type[seedling_indices]    <-  "seedling"
  sim$data$alive[seedling_indices]   <-  TRUE
  sim$data$age[seedling_indices]     <-  0
  
    # Seed dispersal process - adds coordinates and mother_id   
    disperseSeed(sim, adult_indices, seedling_indices)
    dispersePollen(sim, adult_indices, seedling_indices)

}


# Determine which plants survive
procSurvival <- function(sim){
    # If it's the first generation, skip
  if(sim$counter$step == 0) next
  
    # If out of bounds and boundary setting is 'unsuitable', plants die
      # Could maybe speed up by only looking at plants that are currently alive
  if(sim$params$boundary_setting == "unsuitable"){
    sim$data$alive[sim$data$pos_x > sim$params$x_max | sim$data$pos_x < 0] <- FALSE
    sim$data$alive[sim$data$pos_y > sim$params$y_max | sim$data$pos_y < 0] <- FALSE 
  }
    
    
    # Make a vector of whether adults that are currently alive will die
  n_adults_alive <-  sum(with(sim$data, alive[type == "adult" & alive == TRUE]),
                         na.rm = TRUE)
  adult_fate  <- sample(c(TRUE, FALSE), n_adults_alive, replace = TRUE,
                        prob = c(sim$params$adult_survival,
                                 1 - sim$params$adult_survival))
  
  sim$data$alive[with(sim$data, which(type == "adult" & alive == TRUE))] <- adult_fate
  
    ## Seedling survival 

    ## If seedling survival is random...
  if(sim$params$seedling_survival == "random"){
    
    n_seedlings_alive <- sum(with(sim$data, alive[type == "seedling" & alive == TRUE]),
                             na.rm = TRUE)
    
    seedling_fate  <- sample(c(TRUE, FALSE), n_seedlings_alive, replace = TRUE , 
                        prob = c(sim$params$seedling_survival_prob, 
                                 1 - sim$params$seedling_survival_prob))
    
  sim$data$alive[with(sim$data, which(type == "seedling" & alive == TRUE))] <- seedling_fate
  }
  
  ## If seedling survival is distance dependent...
  if(sim$params$seedling_survival == "distance"){
    
    # Find coords of alive seedlings and adults
  seedling_coords <- sim$data[with(sim$data, which(type == "seedling" & alive == TRUE)),
                              c("pos_x", "pos_y")]
  adult_coords <- sim$data[with(sim$data, which(type == "adult" & alive == TRUE)), 
                           c("pos_x", "pos_y")]
  
  # Find distance from each seedling to nearest adult
  distances <- as.matrix(pdist::pdist(seedling_coords, adult_coords)) # calculate distances
  
  dist_to_nearest_adult <- apply(distances, 1, min)

  # Use distance to estimate probability of survival in a logistic form?
  
  prob_of_survival_dist <- plogis(dist_to_nearest_adult,
                             location = sim$params$seedling_survival_dist_location,
                             scale = sim$params$seedling_survival_dist_scale)
  # 1 is alive, 0 is dead
  seedling_fate  <- rbinom(length(prob_of_survival_dist), 1, prob_of_survival_dist)
  seedling_fate <- seedling_fate == 1
  
  sim$data$alive[with(sim$data, which(type == "seedling" & alive == TRUE))] <- seedling_fate
  
  }
  
  
}

# Increase age by one step
increaseAge <- function(sim){
  alive <- which(sim$data$alive == TRUE)
  sim$data$age[alive] <- sim$data$age[alive] + 1
}

# Transition to higher age class
transitionType <- function(sim){
  sim$data$type[sim$data$type == "seedling" & 
                  sim$data$age >= sim$params$age_at_adult] <- "adult"
}







