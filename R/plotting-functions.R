
  # Plot simulation world
plotSim <- function(sim, step = sim$counter$step, alpha_value = 1){
  
    # Initialize plot
   plot(0, type = 'n', asp = 1,
        pch = 19, xlab = "", ylab = "", las = 1, 
        xlim = c(0, sim$params$x_max), ylim = c(0, sim$params$y_max),
        )
   
    # Draw border of study area
   rect(xleft = 0, ybottom = 0, xright = sim$params$x_max, 
        ytop = sim$params$y_max, lty = 2)
   
    # Subset alive seedlings and adults
   adults_alive <- sim$data[sim$registry[[step]][[1]], ]
   seedlings_alive <- sim$data[sim$registry[[step]][[2]], ]
   
    # Plot seedlings, 
    points(seedlings_alive$pos_x, seedlings_alive$pos_y, pch = 21, cex = .75,
            bg = scales::alpha(seedlings_alive$color, alpha_value),
            col = "grey4")
   
   # Plot adults
    points(adults_alive$pos_x, adults_alive$pos_y, pch = 22, cex = 1.25, 
           bg = scales::alpha(adults_alive$color, 1),
           col = "grey4")
}

  # Plot summary data
plotSummary <- function(sim){
  
  old_par <- par("mfrow", "mar")
  
    # Set colors
  adult_col <- "grey4"
  seedling_col <- "aquamarine3"
  
  par(mfrow = c(2,3), mar = c(3.1, 4, 2, 1))
    # Plot adult population size
  plot(sim$summary$generation, sim$summary$n_adults_alive, type = "l", lwd = 2,
       las = 1, main = "Adults", ylab = "N", xlab = "", col = adult_col,
       ylim = c(0, max(sim$summary$n_adults_alive)))
  
    # Plot seedling population size
  plot(sim$summary$generation, sim$summary$n_seedlings_alive, type = "l", lwd = 2,
       las = 1, main = "Seedlings", ylab = "", xlab = "Steps", col = seedling_col)
  
  # Plot landscape
  plotSim(sim)
  
    # Plot He
  plot(sim$summary$generation, sim$summary$he_adults_alive, type = "l", lwd = 2,
       las = 1, ylab = "He", xlab = "", main = "He",
       ylim = c(0, 1), col = adult_col)
  lines(sim$summary$generation, sim$summary$he_seedlings_alive, type = "l", lwd = 2,
        col = seedling_col)
  
    # Plot Nae
  plot(sim$summary$generation, sim$summary$nae_adults, type = "l", lwd = 2,
            las = 1, ylim = c(0, sim$params$n_alleles_per_loci), 
            ylab = "Nae", xlab = "", main = "Nae", col = adult_col)
    lines(sim$summary$generation, sim$summary$nae_seedlings, type = "l", lwd = 2,
             col = seedling_col)
  
    # Plot Sp over time
  plot(sim$summary$generation, sim$summary$sp_adults, type = "l", lwd = 2,
       las = 1, ylab = "Sp", xlab = "", main = "Sp",
       ylim = c(0, 0.4), col = adult_col)
  lines(sim$summary$generation, sim$summary$sp_seedlings, type = "l", lwd = 2,
        col = seedling_col)  
  
  
    # Reset par settings
  on.exit(par(old_par))
}



# Plot dispersal kernels
plotKernels <- function(sim){
  
  dat <- data.frame(seed_obs = calcSeedDispDistance(sim),
                    pollen_obs = calcPollenDispDistance(sim),
                    pollen_eff_obs = calcEffectivePollenDispDistance(sim)
                    )
  dat$seed_kernel <- rweibull(dim(dat)[1], shape = sim$params$seed_kernel_shape,
                              scale = sim$params$seed_kernel_scale)
  
  dat$pollen_kernel <- rweibull(dim(dat)[1], shape = sim$params$pollen_kernel_shape,
                                scale = sim$params$pollen_kernel_scale)

    # Save kernel density estimates
    seed_obs <- density(dat$seed_obs)
    seed_kernel <- density(dat$seed_kernel)
    pollen_obs <- density(dat$pollen_obs)
    pollen_eff_obs <- density(dat$pollen_eff_obs)
    pollen_kernel <- density(dat$pollen_kernel)
    
    # Find max y values for both pollen and seed for plotting
    max_seed_y <- max(seed_obs$y, seed_kernel$y)
    max_pollen_y <- max(pollen_obs$y, pollen_eff_obs$y, pollen_kernel$y)
    
    
    # Set colors
    seed_obs_col <- "green"
    seed_kernel_col <- "purple"
    pollen_obs_col <- "purple"
    pollen_eff_obs_col <- "green"
    pollen_kernel_col <- "red"
    alpha_value <- .25
    
    # Set par settings
    old_par <- par("mfrow", "mar")
    par(mfrow = c(2,1), mar = c(3.1, 4, 2, 1))
  
    # Seed dispersal kernel 
  plot(seed_obs, xlab = "Distance (m)", las = 1, main = "Seed dispersal",
       ylim = c(0, max_seed_y * 1.1))
  polygon(seed_obs, col = scales::alpha(seed_obs_col, alpha_value))
    
  lines(seed_kernel)
  polygon(seed_kernel, col = scales::alpha(seed_kernel_col, alpha_value))
    
  legend("topright", c("Observed", "Kernel"), 
         fill = c(scales::alpha(seed_obs_col, alpha_value),                                                       scales::alpha(seed_kernel_col, alpha_value)))
  
    # Pollen dispersal kernel
  plot(pollen_obs, xlab = "Distance (m)", las = 1, main = "Pollen dispersal",
       ylim = c(0, max_pollen_y * 1.1))
  polygon(pollen_obs, col = scales::alpha(pollen_obs_col, alpha_value))
    
  lines(pollen_eff_obs)
  polygon(pollen_eff_obs, col = scales::alpha(pollen_eff_obs_col, alpha_value))
    
  lines(pollen_kernel)
  polygon(pollen_kernel, col = scales::alpha(pollen_kernel_col, alpha_value))
  legend("topright", c("Observed", "Effective", "Kernel"), 
            fill = c(scales::alpha(pollen_obs_col, alpha_value), 
                     scales::alpha(pollen_eff_obs_col, alpha_value),
                     scales::alpha(pollen_kernel_col, alpha_value)))
  
  on.exit(par(old_par))
}



plotPollenTrace <- function(sim, step = sim$counter$step, alpha_value = 0.25){
  
  plotSim(sim, step = step, alpha_value = alpha_value)
  
  
  # Plot arrows connecting offspring to parents

      # Choose a random adult that has successfully pollinated
    id_father <- sample(na.omit(unique(sim$data$id_father)), 1)
    id_seedlings <- sim$data$id[sim$data$id_father %in% id_father]
    id_mother <- sim$data$id_mother[id_seedlings]
    
      # Arrow from father to mother
     arrows(x0 = sim$data$pos_x[id_father],
            y0 = sim$data$pos_y[id_father],
            x1 = sim$data$pos_x[id_mother],
            y1 = sim$data$pos_y[id_mother],
            lwd = 1.5,
            angle = 30,
            length = .1,
            col = sim$data$color[id_father])

        # Arrow from mother to seedling
     arrows(x0 = sim$data$pos_x[id_mother],
            y0 = sim$data$pos_y[id_mother],
            x1 = sim$data$pos_x[id_seedlings],
            y1 = sim$data$pos_y[id_seedlings],
            lwd = 1.5,
            angle = 30,
            length = .1,
            lty = 2,
            col = sim$data$color[id_father])
}



