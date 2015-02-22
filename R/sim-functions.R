

initDataFrame <- function(){
  
  ## Creata blank data frame that will store all information about individuals
  ## Also had date of death and date of birth and loci names in python version
  col_names <- c("id", "type", "alive", "age", "color", "pos_x", "pos_y",
                 "id_mother", "id_father")
  df <- data.frame(matrix(ncol = length(col_names), nrow = 10000))
  colnames(df) <- col_names
  return(df)
}

initAdults <- function(params){
  
  data <- initDataFrame()
  n_adults <- params$n_adults_init
  
  # Ids always start at 1 and go up for first adults
  data$id[1:n_adults] <- 1:n_adults
  data$type[1:n_adults] <- "adult"
  data$alive[1:n_adults] <- TRUE
  data$age[1:n_adults] <- params$age_at_adult
  data$color[1:n_adults] <- sample(colors(), n_adults)
  
  # Placement of Adults - currently random
  data$pos_x[1:n_adults] <- runif(n = n_adults, min = 0, max = params$x_max)
  data$pos_y[1:n_adults] <- runif(n = n_adults, min = 0, max = params$y_max)  
  
  # Update plant counter
  counter <<- initCounter()
  counter$plant <<- n_adults
  
  return(data)
}

# Counter keeps track of things like what generation we are on, 
# how many plants have been created
initCounter <- function(){
  return(list(plant = 0, step = 1))
}




## Move forward x generations
stepForward <- function(data, counter, params, steps = 1){
  while(steps > 0){
    
    counter$step <- counter$step + 1
    data <- survival(data = data, params)
    data <- reproduce(data = data, counter = counter, params = params)
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



### Generate summary of data
generateSummary <- function(registry){
  
  ## returns dataframe
  
  # Get pop sizes by look at lenght of registry for each class
  pop_sizes <- unlist(lapply(registry, FUN = function(x) lapply(x, length)))
  
  # Subsetting named vector of pop sizes
  adults_alive <- pop_sizes[which(names(pop_sizes) == "adults_alive")]
  seedlings_alive <- pop_sizes[which(names(pop_sizes) == "seedlings_alive")]
  
  summary_df <- data.frame(generation = 1:counter$step, 
                           n_adults_alive = adults_alive,
                           n_seedlings_alive = seedlings_alive)
  plotSummary(summary_df)
  return(summary_df)
}




