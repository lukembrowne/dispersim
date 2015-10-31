
# Create a new simulation object and initialize it
initSim <- function(params){
  
  #.checkParams(params)
  counter <- initCounter()
  sp <- initSpecies(params)
  x <- runif(n = params$n_ad, min = 0, max = params$x_max) # Locations..
  y <- runif(n = params$n_ad, min = 0, max = params$y_max)
  gen_data <- initGenData(params)
  #initSummary(sim)
  cat("Simulation successfully initialized! \n\n")
  return(list(counter = counter, gen_data = gen_data, sp = sp, x = x, y = y))
}

## Creata blank data frame that will store all information about individuals
initDataFrame <- function(params){
  col_names <- params$loci_names
  df <- data.frame(matrix(ncol = length(col_names), 
                          nrow = params$n_ad))
  colnames(df) <- col_names
  gen_data <- df
  cat("Data frame initialized...\n")
  return(gen_data)
}


# Counter keeps track of things like what generation we are on, 
# how many plants have been created
initCounter <- function(){
  counter <- list(step = 1,
                  seed_immi = 0,
                  pollen_immi = 0,
                  selfing = 0)
  cat("Counter initialized... \n")
  return(counter)
}


# Place initial adults and assign them genotypes
initSpecies <- function(params){
  
  # Assign species ID
  species <- rep(1:params$n_sp_init, times = params$n_ad/params$n_sp_init)
  
  if(params$n_ad %% params$n_sp_init != 0) {
    stop("Number of initial species must be a factor of total number of adults.. \n")
  }
  
  return(species)
}
  

# Initialize summary
initSummary <- function(sim){
  max_gen = 50
  sim$summary <- data.frame(generation = rep(NA, max_gen), 
                           n_adults_alive = rep(NA, max_gen),
                           n_seedlings_alive = rep(NA, max_gen),
                           he_adults_alive = rep(NA, max_gen),
                           he_seedlings_alive = rep(NA, max_gen),
                           sp_adults = rep(NA, max_gen),
                           sp_seedlings = rep(NA, max_gen),
                           nae_adults = rep(NA, max_gen),
                           nae_seedlings = rep(NA, max_gen)
  )
  cat("Summary initialized... \n")
}


## Move forward x generations
runSim <- function(params){
  
  ## Save out frequently used sim data to save time calling sim$x everywhere
  # Takes ~ 500 nanoseconds to call params$...
  # Takes ~ 40 nanoseconds to call just ...
  n_ad <- params$n_ad
  dispersal_type <- params$dispersal_type
  boundary <- params$boundary
  seed_immi_rate <- params$seed_immi_rate
  pollen_immi_rate <- params$pollen_immi_rate
  x_max <- params$x_max
  y_max <- params$y_max
  seed_disp_dist <- params$seed_disp_dist
  pollen_disp_dist <- params$pollen_disp_dist
  n_alleles_per_loci <- params$n_alleles_per_loci
  selfing_rate <- params$selfing_rate
  n_loci <- params$n_loci
  n_alleles_per_loci <- params$n_alleles_per_loci
  recruit_thresh <- params$recruit_thresh
  neighborhood_size <- params$neighborhood_size
  dist_beta <- params$dist_beta
  gen_beta <- params$gen_beta
  breed_mode <- params$breed_mode
  save_status <- params$save_status
  save_status_freq <- params$save_status_freq
  plotting <- params$plotting
  steps <- params$steps
  min_pop <- params$min_pop
  replicates <- params$replicates
  temp_data_prefix <- params$temp_data_prefix
  
  ## Initialize matrices that save multiple replicates
  sp_rep <- matrix(NA, nrow = replicates, ncol = n_ad)
  x_rep <- matrix(NA, nrow = replicates, ncol = n_ad)
  y_rep <- matrix(NA, nrow = replicates, ncol = n_ad)
  gen_data_rep <- list()
  
  
  ## Initialize matrices that save status over time
  he_over_time <- matrix(NA, nrow = replicates, ncol = save_status_freq + 1)
  sp_over_time <- matrix(NA, nrow = replicates, ncol = save_status_freq + 1)
  
  
for(replicate in 1:replicates){  
  
  ## Initialize simulation
  sim <- initSim(params = params)
  
  over_time = 1 # Use for indexing over_time matrices
  
  # Save to speed up when accessing data later
  x <- sim$x
  y <- sim$y
  sp <- sim$sp
  gen_data <- sim$gen_data
  counter <- sim$counter
  
  cat("--------- Starting Replicate:", replicate, "--------- \n")
  
  for(step in 1:steps){
    
    # Save on initial step and end step as well
    if(step == 1 | (step %% (steps/save_status_freq) == 0)){
      cat("On generation: ", step/n_ad, "... Step:", step, "... \n", sep = "")
      
      if(plotting == TRUE) {
        plotLandscape(x = x, y = y, sp = sp, 
                      x_max = x_max, y_max = y_max)
      }
      
      if(save_status == TRUE){
        
        sp_abund <- table(sp)
        
        # Species richness
         sp_over_time[replicate, over_time] <- length(sp_abund)
         
         # Calculate expected heterozygosity by species and average it
         he_sp = NULL
         he_ind <- 1
         species_list <- names(sp_abund[sp_abund > min_pop])
         for(species in species_list){
           he_sp[he_ind] <- calcHeAvg(gen_data[sp == species, ], n_loci, n_alleles_per_loci)
           he_ind <- he_ind + 1
         }
         
         he_over_time[replicate, over_time] <- mean(he_sp)
        
         over_time = over_time + 1
      } # End save status
    } 
    
      counter$step <- counter$step + 1
      
      flag_survival = 0  # Set a flag that will get flipped when an offspring survives
      
      # Randomly chose an adult to die, save that index 
      dead_index <- sample(1:n_ad, 1)
      
      # Set important values to NA to ensure it doesn't affect future calculations
      sp[dead_index] <- NA
      x[dead_index] <- NA
      y[dead_index] <- NA
      gen_data[dead_index, ] <- NA
    
      ## Calculate probability of immigration
      seed_immi_prob <- runif(1, min = 0, max = 1)
      
      # Keep track of how many tries we are doing
        try_attempt = 0 
      
      ## Reiterate dispersal and such until survival happens
      while(flag_survival == 0){ 
        
        try_attempt = try_attempt + 1 ## Increment number of tries
        
        if(try_attempt %% 50 == 0 & try_attempt > 0) cat("On attempt #:", try_attempt, "\n")
        
        ## If no immigration... 
        # Divide by n_ad to standardize for generation time
        if(seed_immi_prob > (seed_immi_rate / n_ad)){ 
          
          ## Choose parent to reproduce
          seed_source_index <- sample(1:n_ad, size = 1)
          
          ## If seed source is dead.. try again
          while(seed_source_index == dead_index){
            seed_source_index <- sample(1:n_ad, size = 1)
          }
          
          #### BEGIN DISPERSAL PROCESS
          
          ## If global dispersal, offspring location is random
          if(dispersal_type == "global"){
            offspr_x <- runif(1, min = 0, max = x_max)
            offspr_y <- runif(1, min = 0, max = y_max)
            
          } else if(dispersal_type == "local"){
            
            ## Choose dispersal distance
            # Here is a exponential distribution
            distance <- rexp(1, rate = 1/seed_disp_dist)
            angle <- runif(1, min = 0, max = 2*pi)   
            
            ## Assign new XY coordinates for Offspring
            
            offspr_x <- x[seed_source_index] + cos(angle)*distance  
            offspr_y <- y[seed_source_index] + sin(angle)*distance
          }
      
          if(boundary == "torus"){
            # Eliminate edge effects by simulating a torus
            # Code from Sedio and Ostline 2013 Ecol. Letter
            while(offspr_x < 0) {offspr_x = x_max - offspr_x} 
            while(offspr_x > plot_max) {offspr_x = offspr_x - x_max}
            while(offspr_y < 0) {offspr_y = y_max - offspr_y}
            while(offspr_y > plot_max) {offspr_y = offspr_y - y_max}
          } 
          
          if(boundary == "edge"){
            # Choose a new adult to reproduce if offspring lands out of bounds
            if(offspr_x < 0 | offspr_x > x_max | offspr_y < 0 | offspr_y > y_max){
              next
            }
            
          }
          
          # Check to make sure offspring landed in bounds
          # if(offspr_x < 0 | offspr_y < 0 | offspr_x > plot_max | offspr_y > plot_max){
          #   stop("Offspring landed out of bounds...\n")
          # }
          
          ## Assign species ID (inherits from parent)
          offspr_sp<- sp[seed_source_index]
          
          ## Assign genotype to offspring
          ## Calculate probability of pollen immigration 
          pollen_immi_prob <- runif(1, min = 0, max = 1)
          
          ## If pollen immigration - assign new random genotype
          if(pollen_immi_prob < (pollen_immi_rate / n_ad)){
            
            offspr_gen <- sample(1:n_alleles_per_loci, size = ncol(gen_data), replace = TRUE)
            
          } else
            
            ## Clonal reproduction - directly inherits genotype from mother
            if(breed_mode == "clonal"){
              offspr_gen <- gen_data[seed_source_index, ]
              
            } else 
              
              ## Outbreeding - genotype is the average of mother and father
              ## Includes some chance of self-fertizliation - 1/n_ad
              
              if(breed_mode == "outbreeding") {
                
                self_prob <- runif(1, min = 0, max = 1)
                ## Check for selfing...
                if(self_prob < selfing_rate){
                  
                  offspr_gen <- gen_data[seed_source_index, ] # Set genotype the same as mother
                  
                } else {
                  
                  ## Find indices of the same species
                  potential_father_indices <- which(sp == offspr_sp)
                  
                  ## If potential father is dead... remove from list
                  potential_father_indices <- potential_father_indices[potential_father_indices 
                                                                       != dead_index]
                  
                  if(dispersal_type == "global"){
                    
                    if(length(potential_father_indices) == 1 ){
                      father_index <- potential_father_indices
                    } else {
                      father_index = sample(potential_father_indices, size = 1) 
                    }
                  }
                  
                  # If there's only one potential father and self, otherwise sample from group
                  if(length(potential_father_indices) == 1 ){
                    
                    father_index <- potential_father_indices
                    
                  } else if(dispersal_type == "local"){
                    
                    ### Calculate probability of pollen dispersal using exponential function 
                    # and distance to potential dads
                    if(boundary == "torus"){
                      dist_to_fathers <- .distTorus(x1 = x[seed_source_index], 
                                                   x2 = x[potential_father_indices],
                                                   y1 = y[seed_source_index], 
                                                   y2 = y[potential_father_indices],
                                                   xmax = x_max, ymax = y_max)
                    }
                    
                    if(boundary == "edge"){
                      dist_to_fathers <- .calcDist(x[seed_source_index],
                                                   x[potential_father_indices],
                                                   y[seed_source_index],
                                                   y[potential_father_indices] )
                    }
                    
                    ## Calculate probability of pollination based on pollen dispersal kernel
                    probs_of_pollination <- dexp(dist_to_fathers, rate = 1 / pollen_disp_dist)
                    
                    ## To Standardize to 1 to offset pollen limitation?
                    ## At least a father is chosen every time
                    ## Will tend to increase realized pollen dispersal distances
                    father_index = sample(x = potential_father_indices, size = 1,
                                          prob = probs_of_pollination)
                  }
                  
                  ## Error check - make sure mother and father are the same species
                  if(sp[seed_source_index] != sp[father_index]) {
                    stop("Mom and dad not same species...\n")
                  }
                  
                  ## At this point - if selfing is only option, jump to next while loop 
                  # Because we've already passed selfing filter based on rate
                  if(father_index == seed_source_index){
                    next
                  }
                  
                  # Breed mother and father together to produce offspring genotype
                  offspr_gen <- chooseOffsprGen(mom_gen = gen_data[seed_source_index,],
                                  dad_gen = gen_data[father_index, ],
                                  n_loci = n_loci)
                  
                }  # End finding father ELSE statement
              } ## End outbreeding 
          
        } else 
          ## If immigration from regional pool happens
          if(seed_immi_prob <= (seed_immi_rate / n_ad)){  
            
            ## Assign species ID
            ## Chose random species from regional pool
            offspr_sp <- sample(1:params$n_sp_pool, size = 1) 
            
            ## Assign random spatial coordinates
            offspr_x <- runif(1, min = 0, max = x_max)
            offspr_y <- runif(1, min = 0, max = y_max)
            
            ## Assign random genotype
            offspr_gen <- sample(1:n_alleles_per_loci, size = ncol(gen_data), replace = TRUE)
            
          } ## End immigration IF
        
        
        ## Calculate population level allele frequency for just species in question
        # Initialize matrix
        # Could probably speed up with RCPP
    
        if(gen_beta != 0){
          
          same_sp_indices <- which(sp == offspr_sp)
          gen_data_sub <- gen_data[same_sp_indices, , drop = FALSE]

          ref_al_freq <- calcAlleleFreqCpp(gen_data = gen_data_sub,
                                            n_loci = n_loci,
                                            n_alleles_per_loci = n_alleles_per_loci)
        } else {
          ref_al_freq = NULL
        }
   
        ## Check to see if offspring will live or die
        prob_survival <- procSurvival(offspr_x = offspr_x, offspr_y = offspr_y,
                                           offspr_sp = offspr_sp, 
                                           offspr_gen = offspr_gen,
                                           x = x, y = y, 
                                           sp = sp,
                                           gen_data = gen_data, 
                                           x_max = x_max, y_max = y_max,
                                           recruit_thresh = recruit_thresh,
                                           dead_index = dead_index,
                                           neighborhood_size = neighborhood_size,
                                           dist_beta = dist_beta,
                                           gen_beta = gen_beta,
                                           boundary = boundary,
                                           n_loci = n_loci,
                                           n_alleles_per_loci = n_alleles_per_loci,
                                           ref_al_freq = ref_al_freq)
        
        ## Choose whether offspring dies or not
        
        ## If offspring survives...
        if(rbinom(1, size = 1, prob = prob_survival) == 1) { 
          
          ## Assign information to dead adult index
          x[dead_index] <- offspr_x
          y[dead_index] <- offspr_y
          sp[dead_index] <- offspr_sp
          gen_data[dead_index, ] <- offspr_gen
          
          flag_survival = 1 ## Get out of while loop
          
          # Increment counters
          if(seed_immi_prob <= (seed_immi_rate / n_ad)){
            counter$seed_immi = counter$seed_immi + 1
            counter$pollen_immi = counter$pollen_immi + 1 # Immigrant seed counts as pollen imm
          }
          
          if(pollen_immi_prob < (pollen_immi_rate / n_ad)){
           counter$pollen_immi = counter$pollen_immi + 1
          }
          
          if(self_prob < selfing_rate){
            counter$selfing = counter$selfing + 1
          }
          
        } # End offspring survival IF
      } ## End flag_survival while loop
    } # End step loop
    
    
    cat("Total # of selfing events:", counter$selfing, "...",
        round(counter$selfing/steps, 4), "% per successful breeding event \n")
    cat("Total # of seed immigrants:", counter$seed_immi, "...",
        round(counter$seed_immi/steps*n_ad, 4), "% per generation \n")  
    cat("Total # of pollen immigrants:", counter$pollen_immi, "...",
        round(counter$pollen_immi/steps*n_ad, 4), "% per generation \n")
    
    cat("--------- Ending Replicate:", replicate, "--------- \n")
    cat("\n\n")
    
    
    ### Save output after final generation to summary matrix
    
    sp_rep[replicate, ] <- sp
    x_rep[replicate, ] <- x
    y_rep[replicate, ] <- y
    gen_data_rep[[replicate]] <- gen_data
    
    ## 
    
    temp_list <- list(params = params, 
                      steps = steps,
                      sp_rep = sp_rep,
                      x_rep = x_rep,
                      y_rep = y_rep,
                      gen_data_rep = gen_data_rep,
                      sp_over_time = sp_over_time,
                      he_over_time = he_over_time)
    save(temp_list, 
         file = paste("temp_", temp_data_prefix,"_rep_", replicate, ".Rdata", sep = ""))
    
    

} # End replicate loop    
  
  return(list(params = params, 
              steps = steps,
              sp_rep = sp_rep,
              x_rep = x_rep,
              y_rep = y_rep,
              gen_data_rep = gen_data_rep,
              sp_over_time = sp_over_time,
              he_over_time = he_over_time))
}   # End function     



