

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


