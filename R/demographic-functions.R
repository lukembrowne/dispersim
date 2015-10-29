

# Returns probability of survive for a given offspring
procSurvival <- function(offspr_x, offspr_y, 
                         offspr_gen, offspr_sp,
                         x, y, gen_data, sp,
                         x_max, y_max, 
                         recruit_thresh, 
                         dead_index, 
                         neighborhood_size,
                         dist_beta, gen_beta, 
                         boundary,
                         n_loci, n_alleles_per_loci,
                         ref_al_freq){
  
  # Check if over recruitment threshold
  if(boundary == "torus"){
    distances_all_adults <- .distTorus(x1 = offspr_x, x2 = x, 
                                      y1 = offspr_y, y2 = y,
                                      xmax = plot_max, ymax = plot_max)
  }
  
  if(boundary == "edge"){
    distances_all_adults <- .calcDist(offspr_x, x,
                                      offspr_y, y)
  }
  
  # If below threshold, automatic death
  if(any(distances_all_adults <= recruit_thresh, na.rm = TRUE)){ 
    return(surv_prob = 0)
  }
  
  # Skip costly calculations if distance and genetics don't matter..
  if(dist_beta != 0 | gen_beta != 0){
    
    ### Find indices of the same species
    same_species_indices <- which(sp == offspr_sp)
    
    # If no other adults of that species...
    if(length(same_species_indices) == 0){
      n_neighborhood = 0 
      genetic_effect = 0 
      
    } else {
      
      ## Remove dead tree from same_species_indices
      same_species_indices <- same_species_indices[same_species_indices != dead_index]
      
      # Calculate distances to same species
      distances_to_sp <- distances_all_adults[same_species_indices]
      
      # Find neighbors within neighborhood
      in_neighborhood_indices <- same_species_indices[distances_to_sp < neighborhood_size]
      n_neighborhood <- length(in_neighborhood_indices)
      
      if(n_neighborhood == 0){
        
        genetic_effect = 0 
        
      } else if(gen_beta != 0){
       
        genetic_effect <- calcFij(offspr_gen, gen_data[in_neighborhood_indices, , drop = FALSE],
                                      ref_al_freq,
                                      n_loci = n_loci, 
                                      n_alleles_per_loci = n_alleles_per_loci,
                                      n_gene_copies = length(same_species_indices)*2)
        
        genetic_effect <- mean(genetic_effect)
        
        if(genetic_effect < 0) genetic_effect = 0 # Avoid negative relatedness

      } else if(gen_beta == 0) {
        
        genetic_effect = 0
        
      }
    }
  }
  
  if(dist_beta == 0 & gen_beta == 0){
    n_neighborhood = 0
    genetic_effect = 0
  }
  
  ## Calculate total survival probability
  # We take the mean instead of sum of genetic_effect makes it density independent
  surv_prob <- 1 / (1 + n_neighborhood * dist_beta + 
                      genetic_effect * gen_beta)
  
  return(surv_prob)
}
  
  







