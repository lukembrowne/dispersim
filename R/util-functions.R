# Simple function to check if number is a positive integer
isPositiveInteger <- function(number){
  (number > 0) & (number %% 1 == 0)
}


# Function to check that parameters are within expected ranges
checkParams <- function(params){
  
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
  if(params$boundary_setting != "unsuitable"){
    cat("Unrecognized boundary setting...\n")
    n_errors <- n_errors + 1
  }
  if(!isPositiveInteger(params$n_adults_init)){
    cat("n_adults_init must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(!isPositiveInteger(params$n_loci)){
    cat("n_loci must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(!isPositiveInteger(params$n_alleles_per_loci)){
    cat("n_alleles_per_loci must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(!is.character(params$loci_names)){
    cat("loci_names must be in character format...\n")
    n_errors <- n_errors + 1
  }
  if(!isPositiveInteger(params$age_at_adult)){
    cat("age_at_adult must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(!isPositiveInteger(params$crop_size)){
    cat("crop_size must be positive integer...\n")
    n_errors <- n_errors + 1
  }
  if(params$adult_survival < 0 | params$adult_survival > 1){
    cat("adult_survival must be between 0 and 1...\n")
    n_errors <- n_errors + 1
  }
  if(params$seedling_survival < 0 | params$seedling_survival > 1){
    cat("seedling_survival must be between 0 and 1...\n")
    n_errors <- n_errors + 1
  }
  
  if(n_errors > 0 ){
    stop("Errors detected in parameter list... please fix and try again! \n")
  }
  
  cat("No errors detected in parameter list... Everything looks great.\n")
}
