
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


## Calculate unbiased gene diversity for one locus
## Formula from Nei 1987 that works well with both haploid and diploid data
## Returns mean He across all loci
calcHe <- function(al_freq, n_ind){
    
    # Formula for unbiased gene diversity from Nei 1987
  he <- apply(al_freq, 1, function(x) ((2 * n_ind) / (2 * n_ind - 1)) * (1 - sum(x^2)))

    return(sum(he)/length(he))
}



addNaeToSummary <- function(sim){
  for(step in 1:sim$counter$step){
    sim$summary$nae_adults[step] <- extractDivParam(sim$spagedi_data_adults[[step]], "nae")[1, 1]
    
    if(step == 1) next # Skip first step for seedlings
    sim$summary$nae_seedlings[step] <- extractDivParam(sim$spagedi_data_seedlings[[step]], "nae")[1, 1]
  }
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

  
