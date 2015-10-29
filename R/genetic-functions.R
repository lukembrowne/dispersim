
## Initialize loci for beginning adults
initGenData <- function(params){
  
  col_names <- params$loci_names
  gen_data <- data.frame(matrix(ncol = length(col_names), 
                          nrow = params$n_ad))
  colnames(gen_data) <- col_names

  # Fill in loci
  gen_data <- apply(gen_data, 2, 
                   function(x) {
                       x <- sample(1:params$n_alleles_per_loci, length(x), replace = TRUE)
                        })
  cat("Loci for Adults initialized...\n")
  
  return(gen_data)
}


# Generate Loci names using L1_1 and L1_2 for diploid tissue
makeLociNames <- function(n_loci){ 
  paste('L', as.character(rep(1:n_loci, each = 2)), c("_1", "_2"), sep = "") 
}



# Choose one allele from each diploid loci to pass on to offspring
# Returns new offspring genotype
chooseOffsprGen <- function(mom_gen, dad_gen, n_loci){
  
  base_indices <- seq(from = 1, to = n_loci * 2, by = 2)
  mom_gen_indices <- base_indices + sample(c(0,1), n_loci, replace = TRUE)
  dad_gen_indices <- base_indices + sample(c(0,1), n_loci, replace = TRUE)
   
  offspr_gen <- c(rbind(mom_gen[mom_gen_indices], dad_gen[dad_gen_indices]))
   
  return(offspr_gen)

}

## Calculate allele frequencies
# Returns a named numeric vector
calcAlleleFreq <- function(alleles_1, alleles_2 = NULL, n_alleles_per_loci){
  
  ## Initialize blank al frequency table with all possible alleles
  
  freq <- rep(0, n_alleles_per_loci)
  
    # Get rid of NAs in case they are there
  alleles <- c(alleles_1, alleles_2)
  alleles <- na.omit(alleles)
  
  
    # Sum number of instances per allele and divide by total length
  allele_tab <- table(alleles)
  ind_freq <- as.numeric(allele_tab)/ length(alleles)

  freq[as.numeric(names(allele_tab))] <- ind_freq
  
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

## Calculate He across loci, input is just genotypes of samples you want to analyze

calcHeAvg <- function(sim, data_subset){
        # Could speed up by setting this outside the function somewhere
  loci_names_single <- sim$params$loci_names[grepl("_1", sim$params$loci_names)]
  
    # Init holder to hold he values for each locus
  he_holder <- rep(NA, length(loci_names_single))
    # Loop through and calc He for each locus, then average across loci
  i <- 1
  for(locus in loci_names_single){
    dip_gen <- getDiploidGenotype(data_subset, locus_name = locus)
    he_holder[i] <- calcHe(dip_gen)
    i <- i + 1
  }
  return(mean(he_holder))
}


# Function to return full diploid genotype given only the locus name
# Give just the locus name without the underscore
# Returns a vector of alleles
getDiploidGenotype <- function(data_subset, locus_name, sep = "-"){
  
  if(!is.character(locus_name)){
    stop("Locus name must be formatted as string")
  }
  
  if(!locus_name %in% names(data_subset)){
    stop("Name of locus not found in data")
  }
  
  # Find column index for locus, and also next locus that contains 2nd half 
  # of diploid data
  first_col_index <- grep(locus_name, names(data_subset))
  second_col_index <- first_col_index + 1
  
  if(length(first_col_index) > 1){
    stop("More than one locus found with that name.. please be more specific")
  }
  
  alleles <- c(data_subset[, first_col_index], data_subset[, second_col_index])
    
  return(alleles)
}



addNaeToSummary <- function(sim){
  for(step in 1:sim$counter$step){
    sim$summary$nae_adults[step] <- extractDivParam(sim$spagedi_data_adults[[step]], "nae")[1, 1]
    
    if(step == 1) next # Skip first step for seedlings
    sim$summary$nae_seedlings[step] <- extractDivParam(sim$spagedi_data_seedlings[[step]], "nae")[1, 1]
  }
}


## Calc Fij kinship coefficient from Loiselle et al. 1995
# Between offspring and all eligible adults
# Formula taken from Spagedi 1.4 manual
## Returns mean Fij estimate to each adult

calcFij <- function(offspr_gen, gen_data_sub, ref_al_freq_sub, n_loci, n_alleles_per_loci,
                    n_gene_copies){
  
  fij = NULL
  
  loci_seq <- seq(1, n_loci*2, 2) # Used to loop over loci
  
  ## Calculate individual allele frequency for offspring
  offspr_al_freq <- matrix(0, nrow = n_loci, ncol = n_alleles_per_loci)
  col = 1
  for(locus in 1:n_loci){
    offspr_al_freq[locus, ] <- calcAlleleFreq(offspr_gen[col], offspr_gen[col+1], 
                                              n_alleles_per_loci)
    col = col+2
  }
  
  
  ## Loop through adults
  for(ad in 1:nrow(gen_data_sub)){
    
    ## Calculate individual allele frequency for adult
    ad_al_freq <- matrix(0, nrow = n_loci, ncol = n_alleles_per_loci)
    col = 1
    for(locus in 1:n_loci){
      ad_al_freq[locus, ] <- calcAlleleFreq(gen_data_sub[ad, col], 
                                            gen_data_sub[ad, col+1], 
                                                n_alleles_per_loci)
      col = col+2
    }
    
    
    denom = 0 # Initialize denominator and numerator
    numer = 0
    # Loop through loci
    for(locus in 1:n_loci){
      for(allele in 1:n_alleles_per_loci){
        
        ## Numerator calcuations
        numer_temp =  sum((offspr_al_freq[locus, allele] - ref_al_freq_sub[locus, allele]) *
                            (ad_al_freq[locus, allele] - ref_al_freq_sub[locus, allele])) + 
          sum(ref_al_freq_sub[locus, allele]*(1 - ref_al_freq_sub[locus, allele]) / (n_gene_copies - 1))
        numer = numer + numer_temp
        # Denominator calculations
        denom_temp = ref_al_freq_sub[locus, allele] * (1 - ref_al_freq_sub[locus, allele])
        denom = denom_temp + denom
        
      } # End allele loop
    } # End locus loop
    
    fij[ad] = numer / denom ## Average across loci
  
  } # End adult loop
  
   return(fij)
}





# Calculate average Fij across a population
calcFijPopulation <- function(sim, ids){
  
  ref_gen <- sim$data[ids, sim$params$loci_names] # Calculate
  n_gene_copies <- length(ids) * 2  # Gene copies (* 2 for diploid)
  
  
  id_combos <- combn(ids, 2)
  fijs <- rep(NA, dim(id_combos)[2])
  i <- 0
  for(column in 1:dim(id_combos)[2]){
    fijs[column] <- calcFijPairwise(sim, id1 = id_combos[1, column], 
                                    id2 = id_combos[2, column],
                                     ref_gen = ref_gen, n = n_gene_copies)
  }
  
  return(mean(fijs))
}





calcNae <- function(sim, ids){
  # Following formula (16) in Nielsen et al. 2003
  gen_data <- sim$data[ids, sim$params$loci_names]
  naes <- rep(NA, sim$params$n_loci)
  n_samples <- length(ids) * 2 # Number of gene copies (*2 for diploid)
  i <- 1
  
  for(locus in seq(1, length(sim$params$loci_names), 2)){ # Repeat of Fij pairwise 
    ref_freq <- calcAlleleFreq(gen_data[, locus], gen_data[, locus + 1])
    
    nae <- ((n_samples - 1) ^2) / (sum(ref_freq^2) * (n_samples + 1) * (n_samples - 2) + 3 - n_samples)
    
    naes[i] <- nae
    i <- i + 1
  }
  return(mean(naes))
}

  
