
  # Plot simulation world
plotWorld <- function(data, params){
  
    # Initialize plot
   plot(0, type = 'n', asp = 1,
        pch = 19, xlab = "", ylab = "", las = 1, 
        xlim = c(0, params$x_max), ylim = c(0, params$y_max),
        )
   
    # Draw border of study area
   rect(xleft = 0, ybottom = 0, xright = params$x_max, ytop = params$y_max, lty = 2)
   
   alive <- subset(data, data$alive == TRUE)
   
    # Plot adults
   with(alive[alive$type == "adult", ], 
        points(pos_x, pos_y, pch = 22, cex = 1.25, bg = color))
   
    # Plot seedlings
   with(alive[alive$type == "seedling", ], 
        points(pos_x, pos_y, pch = 21, cex = .75, bg = color))
}

  # Plot summary data of population sizes
plotSummary <- function(summary_df, data, params){
  
  old_par <- par("mfrow", "mar")
  
    # Set colors
  adult_col <- "grey4"
  seedling_col <- "red"
  
  par(mfrow = c(2,2), mar = c(3.1, 4, 2, 1))
    # Plot adult population size
  plot(summary_df$generation, summary_df$n_adults_alive, type = "l", lwd = 2,
       las = 1, main = "Adults", ylab = "N", xlab = "", col = adult_col,
       ylim = c(0, max(summary_df$n_adults_alive)))
    # Plot seedling population size
  plot(summary_df$generation, summary_df$n_seedlings_alive, type = "l", lwd = 2,
       las = 1, main = "Seedlings", ylab = "", xlab = "Steps", col = seedling_col)
  
    # Plot He
  plot(summary_df$generation, summary_df$he_adults_alive, type = "l", lwd = 2,
       las = 1, ylab = "He", xlab = "", main = "He",
       ylim = c(0, 1), col = adult_col)
  lines(summary_df$generation, summary_df$he_seedlings_alive, type = "l", lwd = 2,
        col = seedling_col)

  plotWorld(data, params)
    # Reset par settings
  on.exit(par(old_par))
}



# Plot dispersal kernels

plotKernels <- function(data, params){
  
  dat <- data.frame(seed_obs = calcSeedDispDistance(data),
                    pollen_obs = calcPollenDispDistance(data),
                    pollen_eff_obs = calcEffectivePollenDispDistance(data)
                    )
  dat$seed_kernel <- rweibull(dim(dat)[1], shape = params$seed_kernel_shape,
                              scale = params$seed_kernel_scale)
  
  dat$pollen_kernel <- rweibull(dim(dat)[1], shape = params$pollen_kernel_shape,
                                scale = params$pollen_kernel_scale)
  
  dat_melt <- reshape2::melt(dat)
  
  ggplot2::ggplot(dat_melt, ggplot2::aes(x = value, fill = variable)) +
    ggplot2::geom_density(alpha = .2) +
    ggplot2::labs(x = "Distance", y = "") + 
    ggplot2::theme(legend.title=  ggplot2::element_blank())
}



