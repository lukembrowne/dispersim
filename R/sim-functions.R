
# Define sim class - using reference class OOP system
createSimObject <- setRefClass(Class = "Sim",
            fields = list(data = "data.frame",
                          summary = "data.frame",
                          params = "list",
                          registry = "list",
                          counter = "list",
                          spagedi_data_adults = "list",
                          spagedi_data_seedlings = "list"
                          )
            )

# Create a new simulation object and initialize it
initSim <- function(params){
  sim <- createSimObject$new(params = params)
  
  .checkParams(sim$params)
  initDataFrame(sim)
  initCounter(sim)
  initRegistry(sim)
  initAdults(sim)
  initLoci(sim)
  procCensus(sim)
  cat("Simulation successfully inialized! \n")
  return(sim)
}

## Creata blank data frame that will store all information about individuals
initDataFrame <- function(sim){
  ## Also had date of death and date of birth and loci names in python version
  col_names <- c("id", "type", "alive", "age", "color", "pos_x", "pos_y",
                 "id_mother", "id_father", sim$params$loci_names)
  df <- data.frame(matrix(ncol = length(col_names), 
                          nrow = sim$params$expected_total_individuals))
  colnames(df) <- col_names
  sim$data <- df
  return(cat("Data frame initialized...\n"))
}

# Place initial adults and assign them genotypes
initAdults <- function(sim){
  
  n_adults <- sim$params$n_adults_init
  
  # Ids always start at 1 and go up for first adults
  sim$data$id[1:n_adults] <- 1:n_adults
  sim$data$type[1:n_adults] <- "adult"
  sim$data$alive[1:n_adults] <- TRUE
  sim$data$age[1:n_adults] <- sim$params$age_at_adult
  sim$data$color[1:n_adults] <- sample(colorRamps::primary.colors(n_adults), n_adults)
  
  # Placement of Adults - currently random
  sim$data$pos_x[1:n_adults] <- runif(n = n_adults, min = 0, max = sim$params$x_max)
  sim$data$pos_y[1:n_adults] <- runif(n = n_adults, min = 0, max = sim$params$y_max)
  
  # Update plant counter
  sim$counter$plant <- n_adults
  
  cat("Adults initialized... \n")
}

# Counter keeps track of things like what generation we are on, 
# how many plants have been created
initCounter <- function(sim){
  sim$counter <- list(plant = 0, step = 1)
  cat("Counter initialized... \n")
}

# Initialized registry
initRegistry <- function(sim){
  sim$registry <- list()
  cat("Registry initialized... \n")
}


## Move forward x generations
runSim <- function(sim, steps = 1){
  cat("Beginning simulation...", steps, "total steps \n")
  steps_remaining <- steps
  
  while(steps_remaining > 0){
    sim$counter$step <- sim$counter$step + 1
    procSurvival(sim)
    procReproduce(sim)
    increaseAge(sim)
    transitionType(sim)
    updateCounter(sim)
    chooseGenotypesForOffspring(sim)
    procCensus(sim)
    steps_remaining <- steps_remaining - 1
    cat(signif((1 - steps_remaining / steps) * 100, 2), "% complete \n")
  } 
  cat("Updating summary... \n")
  updateSummary(sim)
  cat("--- Simulation complete! ---")
}        


## Saves IDs of all individuals that are alive - to be used with summary function
procCensus <- function(sim){

  sim$registry[[sim$counter$step]] <- list(
    adults_alive = which(sim$data$type == "adult" & sim$data$alive == TRUE),
    seedlings_alive = which(sim$data$type == "seedling" & sim$data$alive == TRUE)
    )
}


# Update plant counter
updateCounter <- function(sim){
  sim$counter$plant <- min(which(is.na(sim$data$id))) - 1
}



### Generate summary of demographic and genetic data
updateSummary <- function(sim){
  
    # Get pop sizes by look at length of registry for each class
  pop_sizes <- unlist(lapply(sim$registry, FUN = function(x) lapply(x, length)))
  
    # Subsetting named vector of pop sizes
  adults_alive <- pop_sizes[which(names(pop_sizes) == "adults_alive")]
  seedlings_alive <- pop_sizes[which(names(pop_sizes) == "seedlings_alive")]
  
    # Calculate He per generation for adults
  he_adults_alive <- rep(NA, length(sim$registry))
      i <- 1
    for(step in 1:length(sim$registry)){
    he_adults_alive[i] <- calcHeAvg(sim = sim, 
                  data_subset = sim$data[unlist(sim$registry[[step]][1]), ])
    i <- i + 1
    }
      
    # Calculate He per generation for seedlings
  he_seedlings_alive <- rep(NA, length(sim$registry))
    i <- 1
    for(step in 1:length(sim$registry)){
      
      if(step == 1){
        i <- i + 1 # Skip first step for seedlings
        next
      } 
      
      he_seedlings_alive[i] <- calcHeAvg(sim = sim,
                                  data_subset = sim$data[unlist(sim$registry[[step]][2]), ])
      i <- i + 1
    }

    # Put together all information into a data frame
  summary_df <- data.frame(generation = 1:sim$counter$step, 
                           n_adults_alive = adults_alive,
                           n_seedlings_alive = seedlings_alive,
                           he_adults_alive = he_adults_alive,
                           he_seedlings_alive = he_seedlings_alive,
                           sp_adults = NA,
                           sp_seedlings = NA
                           )
  sim$summary <- summary_df 
  
    # Plot summary data                  
  plotSummary(sim)
  

}







