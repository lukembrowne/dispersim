

## Move forward x generations
stepForward <- function(data, counter, params, steps = 1){
  while(steps > 0){
    
    counter$step <- counter$step + 1
    data <- survival(data = data, params)
    data <- reproduce(data = data, counter = counter, crop_size = params$crop_size)
    data <- increaseAge(data)
    data <- transitionType(data, params)
    counter <- updateCounter(data = data, counter = counter)
    registry <<- census(data = data, registry = registry, counter = counter)
    steps <- steps - 1 
  } 
  
  counter <<- counter
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

updateCounter <- function(data, counter){
  counter$plant <- min(which(is.na(data$id))) - 1
  return(counter)  
}