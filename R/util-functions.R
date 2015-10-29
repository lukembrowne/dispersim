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

  


#Write data to spagedi output
writeSpagedi <- function(sim, type, step, file_name, dist_int = -5){
    
  if(type == "adult"){
    registry_subset_index <- 1
  } else if(type == "seedling") {
    registry_subset_index <- 2
  } else {
    stop("Type not recognized...\n")
  }
  
  subset_indices<- sim$registry[[step]][[registry_subset_index]]
  
  sink(file_name)
  # Begin first line
  cat(length(subset_indices),     # First - number of samples
      0,                         # Number of category / populations
      2,                         # Number of spatial dimensions
      sim$params$n_loci,        # Number of loci
      1,                         # Number of digits to code alleles
      2,        # Ploidy
      "\n",
      sep = "\t")        # End first line
  
  # Begin second line - Distance intervals
  if(any(dist_int < 0)) cat(dist_int, "\n") # If negative, just put single number
  if(any(dist_int > 0))
  {
    cat(length(dist_int), dist_int, "\n", sep = "\t")
  }
  
  # Begin third line
  cat("id",            # Individual ID
      "pos_x",                    # X coord
      "pos_y",                    # Y coord
      sim$params$loci_names[seq(1, length(sim$params$loci_names), by = 2)],
      "\n",
      sep = "\t")       # End third line
  sink() # End header
  
  # Begin fourth line - DATA!
  df_to_write <- data.frame(id = sim$data$id[subset_indices],
                            pos_x = sim$data$pos_x[subset_indices],
                            pos_y = sim$data$pos_y[subset_indices])
  

  ## Need to paste together allele calls
  
  data_holder <- data.frame(matrix(nrow = length(subset_indices),
                                    ncol = sim$params$n_loci))

  data_holder <- .pasteLoci(sim$data[subset_indices, sim$params$loci_names])  
  
  colnames(data_holder) <- sim$params$loci_names[seq(1, length(sim$params$loci_names), by = 2)]
  
  
  out <- cbind(df_to_write, data_holder)
  
  write.table(out, file_name, append = TRUE, 
              col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  
  cat("END", file = file_name, append = TRUE )
  
}


# Paste allele calls together for output to Spagedi
.pasteLoci <- function(genetic_data){

  genetic_data_pasted <- matrix(NA, nrow = dim(genetic_data)[1], ncol = dim(genetic_data)[2]/2)
  
  j <- 1
  for(col_index in seq(1, ncol(genetic_data), by = 2)){
    pasted_loci <- paste(genetic_data[, col_index], genetic_data[, col_index + 1], sep = "-")
    genetic_data_pasted[, j] <- pasted_loci
    j <- j + 1
  }

  return(genetic_data_pasted)
}











