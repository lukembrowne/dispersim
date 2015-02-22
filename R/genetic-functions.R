
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


# Small function that takes genotypes of mom and dad and creates a new genotype
# for offspring
breed <- function(data, params, seedling_ids){
# Issue here that it can't assign genotypes to seedlings of adults that recruited into the population - need to make the function recursive?
      # Subset out ids of mothers and fathers
    id_mother <- data$id_mother[seedling_ids]
    id_father <- data$id_father[seedling_ids]
    
    # Choose one allele from each maternal and paternal genotype to pass on
  mom_haploid <- t(apply(data[id_mother, params$loci_names], 1, 
                       function(x) '['(x, chooseAllele(params) )))
  dad_haploid <- t(apply(data[id_father, params$loci_names], 1, 
                       function(x) '['(x, chooseAllele(params) )))
  
    # Fill in genotypes for offspring
  fill <- data[seedling_ids, params$loci_names]
  fill[, seq(1, params$n_loci * 2, by = 2)] <- mom_haploid
  fill[, seq(2, params$n_loci * 2, by = 2)] <- dad_haploid
  data[seedling_ids, params$loci_names] <- fill
  
  return(data)
}

# Choose one allele from each diploid loci to pass on to offspring
# These are the indices to use for subseting
chooseAllele <- function(params){
  seq(from = 1, to = params$n_loci * 2, by = 2) +
    sample(c(0,1), params$n_loci/2, replace = TRUE)
}


