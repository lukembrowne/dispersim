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