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

plotSummary <- function(summary_df){
  
  old_par <- par("mfrow", "mar")
  
  par(mfrow = c(2,1), mar = c(3.1, 3, 2, 1))
    # Plot adult population size
  plot(summary_df$generation, summary_df$n_adults_alive, type = "l", lwd = 2,
       las = 1, main = "Adults", ylab = "", xlab = "", 
       ylim = c(0, max(summary_df$n_adults_alive)))
    # Plot seedling population size
  plot(summary_df$generation, summary_df$n_seedlings_alive, type = "l", lwd = 2,
       las = 1, main = "Seedlings", ylab = "", xlab = "Steps")
 
    # Reset par settings
  on.exit(par(old_par))
}