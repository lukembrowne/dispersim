# Simple function to check if number is a positive integer
.isPositiveInteger <- function(number){
  (number > 0) & (number %% 1 == 0)
}


# Calculate distance between two points
.calcDist <- function(x1, x2, y1, y2){
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}


## Calculating distance on a torus from A to B
## Code from Sedio and Ostling 2013 
.distTorus = function(x1, x2, y1, y2, xmax, ymax) 	
{
  xdif = mapply(min, abs(x1-x2), xmax-abs(x1-x2))
  ydif = mapply(min, abs(y1-y2), ymax-abs(y1-y2))
  xydist=sqrt(xdif^2 + ydif^2)
  return(xydist)
} 	# end distfunc


# Function to check that parameters are within expected ranges
.checkParams <- function(params){
  
  n_errors <- 0
  cat("Checking parameter list for errors...\n")
  if(params$x_max <= 0) {
    cat("x_max must be greater than 0...\n")
    n_errors <- n_errors + 1
  }
  if(params$y_max <= 0){
    cat("y_max must be greater than 0...\n")
    n_errors <- n_errors + 1
  }
  if(!(params$boundary %in% c("torus","edge"))){
    cat("Unrecognized boundary setting...\n")
    n_errors <- n_errors + 1
  }
  if(!.isPositiveInteger(params$n_ad)){
    cat("n_ad must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(!.isPositiveInteger(params$n_loci)){
    cat("n_loci must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(!.isPositiveInteger(params$n_alleles_per_loci)){
    cat("n_alleles_per_loci must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(!is.character(params$loci_names)){
    cat("loci_names must be in character format...\n")
    n_errors <- n_errors + 1
  }
  if(!.isPositiveInteger(params$age_at_adult)){
    cat("age_at_adult must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(!.isPositiveInteger(params$crop_size)){
    cat("crop_size must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(params$adult_survival < 0 | params$adult_survival > 1){
    cat("adult_survival must be between 0 and 1...\n")
    n_errors <- n_errors + 1
  }
  if(!(params$seedling_survival %in% c("random", "distance"))){
    cat("seedling_survival must be either 'random' or 'distance'... \n")
    n_errors <- n_errors + 1
  }
  if(params$seedling_survival_prob < 0 | params$seedling_survival_prob > 1){
    cat("seedling_survival must be between 0 and 1...\n")
    n_errors <- n_errors + 1
  }
  
  if(n_errors > 0 ){
    stop("Errors detected in parameter list... please fix and try again! \n")
  }
  
  cat("No errors detected in parameter list... Everything looks great.\n")
}












