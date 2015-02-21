

## Move forward x generations
stepForward <- function(data, counter, params, steps = 1){
  while(steps > 0){
    
    counter$step <- counter$step + 1
    data <- survival(data = data)
    data <- reproduce(data = data, counter = counter, crop_size = params$crop_size)
    registry <<- census(data = data, registry = registry, counter = counter)
    #generate summary
    # print status
    steps <- steps - 1   
  } 
  
  counter$step <<- counter$step
  return(data)
}        


## Saves IDs of all individuals that are alive - to be used with summary function
census <- function(data, registry, counter){
  
  if(!exists("registry")) registry <- list()

  registry[[counter$step]] <- list(
    adults_alive = which(data$type == "adult" & data$alive == TRUE),
    seedlings_alive = which(data$type == "seedling" & data$alive == TRUE))
  
  return(registry)
}