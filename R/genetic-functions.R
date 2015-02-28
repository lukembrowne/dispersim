
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
    
    # Get rid of NAs in case they are there
  alleles <- c(alleles_1, alleles_2)
  alleles <- na.omit(alleles)
  
    # Used in correcting for bias in small sample sizes
    # Will set ploidy_multiplier to 2 if diploid genotypes are passed to function
  ploidy_multiplier <- ifelse(length(alleles_2 > 0), 2, 1)
  
    # Calc number of samples and allele frequencie
  n_samples <- length(alleles_1)
  freq <- calcAlleleFreq(alleles)
  
    # Formula for unbiased gene diversity from Nei 1987
  he <- ((ploidy_multiplier * n_samples) / (ploidy_multiplier * n_samples - 1)) *
    (1 - sum(freq^2))
  
  return(he)
}




