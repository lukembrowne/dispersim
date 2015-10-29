library(shiny)
library(pdist)
library(dispersim)

shinyServer(function(input, output) {
  
  params_reactive <- reactive({
    
    params <- list(
      # Landscape parameters
      x_max = input$x_max,
      y_max = input$y_max,
      boundary_setting = "unsuitable",
      # Expected number of steps and individuals
      expected_total_individuals = 100000,
      # Num of starting adults
      n_ad = input$n_ad,
      # Genetic parameters
      n_loci = input$n_loci,
      n_alleles_per_loci = input$n_alleles_per_loci,
      loci_names = makeLociNames(n_loci = 10),
      # Demographic parameters
      age_at_adult = input$age_at_adult,
      crop_size = input$crop_size,
      adult_survival = input$adult_survival,
      seedling_survival = "random", # Random or distance dependent
      seedling_survival_prob = input$seedling_survival_prob,
      seedling_survival_dist_location = 10,
      seedling_survival_dist_scale = 1,
      # Dispersal parameters
      seed_kernel_scale = input$seed_kernel_scale,
      seed_kernel_shape = input$seed_kernel_shape,
      pollen_kernel_scale = input$pollen_kernel_scale,
      pollen_kernel_shape = input$pollen_kernel_shape
    )
  })
  
  
  ## Density plot of seed dispersal kernel
  output$seed_kernel_plot <- renderPlot({
    
    params <- params_reactive()
    
    sim_seed_kernel <- rweibull(5000, shape = params$seed_kernel_shape,
                                scale = params$seed_kernel_scale)
    
    
    seed_kernel <- density(sim_seed_kernel, na.rm = TRUE)
    
    plot(seed_kernel, xlab = "Distance (m)", las = 1, main = "Seed dispersal",
         xlim = c(0, params$x_max * 1.5))
    polygon(seed_kernel, col = "forestgreen")
    
  
  })
  
  
  output$pollen_kernel_plot <- renderPlot({
    
    params <- params_reactive()
    
    sim_pollen_kernel <- rweibull(5000, shape = params$pollen_kernel_shape,
                                scale = params$pollen_kernel_scale)
    
    
    pollen_kernel <- density(sim_pollen_kernel, na.rm = TRUE)
    
    plot(pollen_kernel, xlab = "Distance (m)", las = 1, main = "pollen dispersal",
         xlim = c(0, params$x_max * 1.5))
    polygon(pollen_kernel, col = "steelblue")
    
    
  })
  
  
  
  
  initialize_simulation_reactive <- reactive({
    
    input$initialize_new_simulation
    params <- params_reactive()
    
    sim <- initSim(params = params)
    
    runSim(sim, steps = input$n_steps)
    
    
    output$sim_landscape_plot <- renderPlot({
      plotSim(sim, color = input$plotting_color)
    })
    
  })

    

  
  output$sim_summary_plot <- renderPlot({
    #input$update_plot
    sim = initialize_simulation_reactive()
  })
  
  
})