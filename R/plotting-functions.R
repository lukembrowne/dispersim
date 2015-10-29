
## General landscape plotting function
plotLandscape <- function(x, y, sp, x_max, y_max){
  
  plot(x, y, asp = 1, pch = 21, type = "n",
       las = 1, xlab = "", ylab = "", xlim = c(0, x_max), ylim = c(0, y_max))
  text(x, y, labels = as.character(sp))
  rect(xleft = 0, ybottom = 0, xright = x_max, ytop = y_max, lty = 2)
}




## Function to summarize output of multiple iterations of the model
summarizeSim <- function(sim_out, quiet = FALSE){
  
  par(mfrow = c(3,2), mar = c(3,5,2,2))
  
  summary <- list()
  
  
  ## Species richness
  
  # Plot number of species over iterations
  if(quiet == FALSE){
    matplot(t(jitter(sim_out$sp_over_time)), type = "b", pch = 19, las = 1, 
            ylab = "Species richness",
            xlab = "Step",
            xaxt = "n")
    axis(side = 1, at = 1:ncol(sim_out$sp_over_time), 
         labels = c(1, (sim_out$params$steps/(ncol(sim_out$sp_over_time)-1) * 
                          1:(ncol(sim_out$sp_over_time)-1)) ))
  }
  
  summary$sp <- apply(sim_out$sp_rep, 1, function(x) length(table(x))) 
  
  if(quiet == FALSE){
    cat("Species richness summary:", mean(summary$sp),"+-", 
        sd(summary$sp)/sqrt(sim_out$params$replicates), "\n")
    boxplot(summary$sp, las = 1, ylab = "Species richness")
    points(jitter(rep(1, length(summary$sp))), summary$sp, pch = 19)
  }
  
  ## Expected heterozygosity over time
  if(quiet == FALSE){
    matplot(t(sim_out$he_over_time), type = "b", pch = 19, las = 1,
            ylab = "He",
            xlab = "Step", xaxt = "n")
    axis(side = 1, at = 1:ncol(sim_out$he_over_time), 
         labels = c(1,(sim_out$params$steps/(ncol(sim_out$sp_over_time)-1) * 
                         1:(ncol(sim_out$sp_over_time)-1))))
  }
    
    
    ## Species abundance plots
    
    for(row in 1:nrow(sim_out$sp_rep)){
      alpha = .55
      if(row == 1){
        barplot(sort(as.numeric(table(sim_out$sp_rep[row, ])), decreasing = TRUE),
                ylab = "Abundance", las = 1, col = rgb(1, 0, 0, alpha),
                xlim = c(0, max(apply(sim_out$sp_rep, 1, function(x) length(table(x))))+1))
      }
      barplot(sort(as.numeric(table(sim_out$sp_rep[row, ])), decreasing = TRUE),
              ylab = "Abundance", las = 1, add = TRUE, 
              col = rgb(runif(1,0,1), runif(1,0,1), runif(1,0,1), alpha),
              main = "Species abundance plot")
    }
    
  on.exit(par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1)))
  return(summary)
  
  }

    


# Plot dispersal kernels
plotKernels <- function(sim, type, step){
  
  if(type == "adult" & step <= sim$params$age_at_adult){
    stop("Cannot calculate dispersal kernel for adults below step ", 
         sim$params$age_at_adult, "...\n")
  }
  
  if(type == "seedling" & step == 1){
    stop("Cannot calculate dispersal kernel for seedlings at step 1...\n")
  }
  
  if(step < 1 | step > sim$counter$step){
    stop("Step is not valid...\n")
  }
  
    # Build dataframe of dispersal distances
  dat <- data.frame(seed_obs = calcDispDistance(sim, type, method = "seed", step),
                    pollen_obs = calcDispDistance(sim, type, method = "pollen", step),
                pollen_eff_obs = calcDispDistance(sim, type, method = "total pollen", step))
  
    # Generate data from set dispersal kernels
  dat$seed_kernel <- rweibull(dim(dat)[1], shape = sim$params$seed_kernel_shape,
                              scale = sim$params$seed_kernel_scale)
  
  dat$pollen_kernel <- rweibull(dim(dat)[1], shape = sim$params$pollen_kernel_shape,
                                scale = sim$params$pollen_kernel_scale)

    # Save kernel density estimates
    seed_obs <- density(dat$seed_obs, na.rm = TRUE)
    seed_kernel <- density(dat$seed_kernel, na.rm = TRUE)
    pollen_obs <- density(dat$pollen_obs, na.rm = TRUE)
    pollen_eff_obs <- density(dat$pollen_eff_obs, na.rm = TRUE)
    pollen_kernel <- density(dat$pollen_kernel, na.rm = TRUE)
    
    # Find max y values for both pollen and seed for plotting
    max_seed_y <- max(seed_obs$y, seed_kernel$y)
    max_pollen_y <- max(pollen_obs$y, pollen_eff_obs$y, pollen_kernel$y)
    max_x <- max(seed_obs$x, seed_kernel$x, pollen_obs$x, pollen_eff_obs$x, pollen_kernel$x)
    
    
    # Set colors
    seed_obs_col <- "green"
    seed_kernel_col <- "grey50"
    pollen_obs_col <- "blue"
    pollen_eff_obs_col <- "red"
    pollen_kernel_col <- "grey50"
    alpha_value <- .25
    
    # Set par settings
    old_par <- par("mfrow", "mar")
    par(mfrow = c(2,1), mar = c(3.1, 4, 2, 1))
  
      # Seed dispersal kernel 
    plot(seed_obs, xlab = "Distance (m)", las = 1, main = "Seed dispersal",
         ylim = c(0, max_seed_y * 1.1), xlim = c(0, max_x * 1.1))
    polygon(seed_obs, col = scales::alpha(seed_obs_col, alpha_value))
      
    lines(seed_kernel)
    polygon(seed_kernel, col = scales::alpha(seed_kernel_col, alpha_value))
      
    legend("topright", c("Observed", "Kernel"), 
           fill = c(scales::alpha(seed_obs_col, alpha_value),                                                       scales::alpha(seed_kernel_col, alpha_value)))
  
    # Pollen dispersal kernel
    plot(pollen_obs, xlab = "Distance (m)", las = 1, main = "Pollen dispersal",
         ylim = c(0, max_pollen_y * 1.1), xlim = c(0, max_x * 1.1))
    polygon(pollen_obs, col = scales::alpha(pollen_obs_col, alpha_value))
      
    lines(pollen_eff_obs)
    polygon(pollen_eff_obs, col = scales::alpha(pollen_eff_obs_col, alpha_value))
      
    lines(pollen_kernel)
    polygon(pollen_kernel, col = scales::alpha(pollen_kernel_col, alpha_value))
    legend("topright", c("Observed", "Total", "Kernel"), 
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

plotSurvivalCurve <- function(sim){
  location <- sim$params$seedling_survival_dist_location
  scale <- sim$params$seedling_survival_dist_scale
  
  plot(1:sim$params$x_max, plogis(1:sim$params$x_max, location = location, scale = scale),
       ylab = "Probability of survival", xlab = "Distance", las = 1, type = "l", lwd = 3,
       ylim = c(0, 1), main = "Distance based survival of seedlings")
}



