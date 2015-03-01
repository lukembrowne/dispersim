
  # Plot simulation world
plotSim <- function(sim){
  
    # Initialize plot
   plot(0, type = 'n', asp = 1,
        pch = 19, xlab = "", ylab = "", las = 1, 
        xlim = c(0, sim$params$x_max), ylim = c(0, sim$params$y_max),
        )
   
    # Draw border of study area
   rect(xleft = 0, ybottom = 0, xright = sim$params$x_max, 
        ytop = sim$params$y_max, lty = 2)
   
   alive <- subset(sim$data, sim$data$alive == TRUE,
                   select = c("type", "color", "pos_x", "pos_y"))
   
    # Plot adults
   with(alive[alive$type == "adult", ], 
        points(pos_x, pos_y, pch = 22, cex = 1.25, bg = color))
   
    # Plot seedlings
   with(alive[alive$type == "seedling", ], 
        points(pos_x, pos_y, pch = 21, cex = .75, bg = color))
}

  # Plot summary data
plotSummary <- function(sim){
  
  old_par <- par("mfrow", "mar")
  
    # Set colors
  adult_col <- "grey4"
  seedling_col <- "red"
  
  par(mfrow = c(2,2), mar = c(3.1, 4, 2, 1))
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



