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







