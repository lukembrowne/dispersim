
## Initialize loci for beginning adults
initLoci <- function(data, params){
    # Select just loci columns to fill in
  fill <- data[1:params$n_adults_init, params$loci_names]
  fill <- apply(fill, 2, 
        function(x) {
          x <- sample(1:params$n_alleles_per_loci, length(x), replace = TRUE)
        })
  data[1:params$n_adults_init, params$loci_names] <- fill
  return(data)
}

# Generate Loci names using L1_1 and L1_2 for diploid tissue
makeLociNames <- function(n_loci){ 
  paste('L', as.character(rep(1:n_loci, each = 2)), c("_1", "_2"), sep = "") 
}


# Function that takes genotypes of mom and dad and creates a new genotype
# for offspring
chooseGenotypesForOffspring <- function(data, params, counter){

      # Find plants that don't have genotypes yet
    where_na <- which(apply(data[1:counter$plant, params$loci_names], 1, anyNA))
    
    if(length(where_na) == 0){
      cat("All plants already have genotypes... exiting chooseGenotypes function")
      return(data)
    }
  
    offspring_ids <- data$id[where_na]
    
      # Subset out ids of mothers and fathers
    id_mother <- data$id_mother[offspring_ids]
    id_father <- data$id_father[offspring_ids]
    
    # Choose one allele from each maternal and paternal genotype to pass on
  mom_haploid <- t(apply(data[id_mother, params$loci_names], 1, 
                       function(x) '['(x, chooseAllele(params) )))
  dad_haploid <- t(apply(data[id_father, params$loci_names], 1, 
                       function(x) '['(x, chooseAllele(params) )))
  
    # Fill in genotypes for offspring
  fill <- data[offspring_ids, params$loci_names]
  fill[, seq(1, params$n_loci * 2, by = 2)] <- mom_haploid
  fill[, seq(2, params$n_loci * 2, by = 2)] <- dad_haploid
  data[offspring_ids, params$loci_names] <- fill
  
    # If there are any plants left that don't have a genotype, breed again!
  if(any(apply(data[1:counter$plant, params$loci_names], 1, anyNA))){
    data <- chooseGenotypesForOffspring(data, params, counter)
  }
  
  return(data)
}

# Choose one allele from each diploid loci to pass on to offspring
# These are the indices to use for subseting
chooseAllele <- function(params){
  seq(from = 1, to = params$n_loci * 2, by = 2) +
    sample(c(0,1), params$n_loci/2, replace = TRUE)
}

## Calculate allele frequencies
# If diploid, give two lists of alleles
calcAlleleFreq <- function(alleles_1, alleles_2 = NULL){
  
    # Get rid of NAs in case they are there
  alleles <- c(alleles_1, alleles_2)
  alleles <- na.omit(alleles)
  
    # Sum number of instances per allele and divide by total length
  freq <- as.numeric(table(alleles) / length(alleles))
  return(freq)
}

## Calculate unbiased gene diversity for one locus
## Formula from Nei 1987 that works well with both haploid and diploid data
calcHe <- function(alleles_1, alleles_2 = NULL){
    
    # Used in correcting for bias in small sample sizes
    # Will set ploidy_multiplier to 2 if diploid genotypes are passed to function
  ploidy_multiplier <- ifelse(length(alleles_2 > 0), 2, 1)
  
    # Calc number of samples and allele frequencie
  n_samples <- length(na.omit(alleles_1))
  freq <- calcAlleleFreq(alleles_1, alleles_2)
  
    # Formula for unbiased gene diversity from Nei 1987
  he <- ((ploidy_multiplier * n_samples) / (ploidy_multiplier * n_samples - 1)) *
    (1 - sum(freq^2))
  
  return(he)
}

## Calculate He across loci

calcHeAvg <- function(data, params){
        # Could speed up by setting this outside the function somewhere
  loci_names_single <- params$loci_names[grepl("_1", params$loci_names)]
  
    # Init holder to hold he values for each locus
  he_holder <- rep(NA, length(loci_names_single))
    # Loop through and calc He for each locus, then average across loci
  i <- 1
  for(locus in loci_names_single){
    dip_gen <- getDiploidGenotype(data, locus_name = locus)
    he_holder[i] <- calcHe(dip_gen)
    i <- i + 1
  }
  return(mean(he_holder))
}


# Function to return full diploid genotype given only the locus name
# Give just the locus name without the underscore
# Returns a vector of alleles
getDiploidGenotype <- function(data, locus_name, mode = "conc", sep = "-"){
  
  if(!is.character(locus_name)){
    stop("Locus name must be formatted as string")
  }
  
  if(!locus_name %in% names(data)){
    stop("Name of locus not found in data")
  }
  
  # Find column index for locus, and also next locus that contains 2nd half 
  # of diploid data
  first_col_index <- grep(locus_name, names(data))
  second_col_index <- first_col_index + 1
  
  if(length(first_col_index) > 1){
    stop("More than one locus found with that name.. please be more specific")
  }
  
  if(mode == "conc"){
     
    alleles <- c(data[, first_col_index], data[, second_col_index])
    
    return(alleles)
  }

  # Add ability to paste the two genotypes together for output in spagedi program
  # Note that this returns a string instead of a numeric
  if(mode == "paste"){
    
    alleles_pasted <- paste(data[, first_col_index], data[, second_col_index],
                            sep = sep)
    return(alleles_pasted)
  }
}




