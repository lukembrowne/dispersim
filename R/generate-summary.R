### Generate summary of data


generateSummary <- function(world){
  
  ## returns dataframe
  
  # Get pop sizes by look at lenght of registry for each class
  pop_sizes <- unlist(lapply(world$registry, FUN = function(x) lapply(x, length)))
  
    # Subsetting named vector of pop sizes
  adults_alive <- pop_sizes[which(names(pop_sizes) == "adults_alive")]
  seedlings_alive <- pop_sizes[which(names(pop_sizes) == "seedlings_alive")]
  
  summary_df <- data.frame(generation = 1:world$step_counter, 
                           n_adults_alive = adults_alive,
                           n_seedlings_alive = seedlings_alive)
  
  return(summary_df)
}

