
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
  seedling_col <- "red"
  
  par(mfrow = c(2,3), mar = c(3.1, 4, 2, 1))
    # Plot adult population size
  plot(sim$summary$generation, sim$summary$n_adults_alive, type = "l", lwd = 2,
       las = 1, main = "Adults", ylab = "N", xlab = "", col = adult_col,
       ylim = c(0, max(sim$summary$n_adults_alive)))
    # Plot seedling population size
  plot(sim$summary$generation, sim$summary$n_seedlings_alive, type = "l", lwd = 2,
       las = 1, main = "Seedlings", ylab = "", xlab = "Steps", col = seedling_col)
  
    # Plot He
  plot(sim$summary$generation, sim$summary$he_adults_alive, type = "l", lwd = 2,
       las = 1, ylab = "He", xlab = "", main = "He",
       ylim = c(0, 1), col = adult_col)
  lines(sim$summary$generation, sim$summary$he_seedlings_alive, type = "l", lwd = 2,
        col = seedling_col)
  
    # Plot Sp over time
  plot(sim$summary$generation, sim$summary$sp_adults, type = "l", lwd = 2,
       las = 1, ylab = "Sp", xlab = "", main = "Sp",
       ylim = c(0, 0.4), col = adult_col)
  lines(sim$summary$generation, sim$summary$sp_seedlings, type = "l", lwd = 2,
        col = seedling_col)  
  
  
    # Plot landscape
  plotSim(sim)
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
  
  dat_melt <- reshape2::melt(dat)
  
  ggplot2::ggplot(dat_melt, ggplot2::aes(x = value, fill = variable)) +
    ggplot2::geom_density(alpha = .2) +
    ggplot2::labs(x = "Distance", y = "") + 
    ggplot2::theme(legend.title=  ggplot2::element_blank())
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



