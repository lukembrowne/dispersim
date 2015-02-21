plotWorld <- function(world){
  
    # Initialize plot
   plot(0, type = 'n', asp = 1,
        pch = 19, xlab = "", ylab = "", las = 1, 
        xlim = c(0, world$x_max), ylim = c(0, world$y_max),
        )
   
    # Draw border of study area
   rect(xleft = 0, ybottom = 0, xright = world$x_max, ytop = world$y_max, lty = 2)
   
    # Plot adults
   with(world$data[world$data$type == "adult", ], 
        points(pos_x, pos_y, pch = 19))
   
    # Plot seedlings
   with(world$data[world$data$type == "seedling", ], 
        points(pos_x, pos_y, pch = 18, cex = .75, col = "green"))
}