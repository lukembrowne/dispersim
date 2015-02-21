
initSim <- function(params){
  
  foo <- list(x_max = x_max,
              y_max = y_max,
              n_adults_init = n_adults_init,
              plant_counter = 0,
              step_counter = 1,
              crop_size = crop_size,
              seed_kernel_scale = seed_kernel_scale,
              seed_kernel_shape = seed_kernel_shape,
              registry = list(),
              data = initDataFrame())
  
  class(foo) <- "world"
  
  foo <- makeAdults(foo)
  foo <- saveIdAlive(foo)
  
  return(foo)
  
}

initDataFrame <- function(){
 
  ## Creata blank data frame that will store all information about individuals
  ## Also had date of death and date of birth and loci names in python version
  col_names <- c("id", "type", "alive", "age", "color", "pos_x", "pos_y",
                 "id_mother", "id_father")
  df <- data.frame(matrix(ncol = length(col_names), nrow = 5000))
  colnames(df) <- col_names
  return(df)
}

makeAdults <- function(world){
  
  n_adults <- world$n_adults_init
  
    # Ids always start at 1 and go up for first adults
  world$data$id[1:n_adults] <- 1:n_adults
  world$data$type[1:n_adults] <- "adult"
  world$data$alive[1:n_adults] <- TRUE
  world$data$age[1:n_adults] <- 0
  world$data$color[1:n_adults] <- sample(colors(), n_adults)
  
    # Placement of Adults - currently random
  world$data$pos_x[1:n_adults] <- runif(n = n_adults, min = 0, max = world$x_max)
  world$data$pos_y[1:n_adults] <- runif(n = n_adults, min = 0, max = world$y_max)  
  
  # Update plant counter
  world$plant_counter <- n_adults
  
  return(world)
}


