library(shiny)
library(dispersim)
library(pdist)

# devtools::install_github("lukembrowne/dispersim") # Run before uploading to Shinyapps.io or will get a cannot find local installation of dispersim error

shinyUI( fluidPage(
 
  fluidRow(
    column(3,
           wellPanel(
             h4("Landscape options"),
      sliderInput("x_max", 
                  "Landscape width:", 
                  value = 100,
                  min = 1, 
                  max = 1000),
      sliderInput("y_max", 
                  "Landscape height:", 
                  value = 100,
                  min = 1, 
                  max = 1000)
         ), # End landscape panel
    
      
      wellPanel(
          h4("Population options"),
          sliderInput("n_ad", 
                      "Number of intial adults:", 
                      value = 50,
                      min = 1, 
                      max = 100),
          sliderInput("crop_size", 
                      "Crop size per adult per step:", 
                      value = 10,
                      min = 1, 
                      max = 100),
          sliderInput("age_at_adult", 
                      "Age for transition to adult:", 
                      value = 5,
                      min = 1, 
                      max = 10)
      ), # End Population options
      
      wellPanel(

         h4("Demographic options"),
         sliderInput("adult_survival", 
                     "Adult survival per step:", 
                     value = 0.9,
                     min = 0, 
                     max = 1),
         
         sliderInput("seedling_survival_prob", 
                     "Seedling survival per step:", 
                     value = 0.4,
                     min = 0, 
                     max = 1)
      ), # End demographic panel 
  
      wellPanel(
         h4("Genetic options"),
         sliderInput("n_loci", 
                     "Number of loci:", 
                     value = 10,
                     min = 1, 
                     max = 100),
         
         sliderInput("n_alleles_per_loci", 
                     "Number of alleles per loci:", 
                     value = 10,
                     min = 1, 
                     max = 30)
      ) # End genetic panel

  ), # End column of options

  column(3,
    
      wellPanel(
        h4("Seed dispersal kernel options"),
        sliderInput("seed_kernel_scale", 
                    "Seed kernel scale:", 
                    value = 10,
                    min = 1, 
                    max = 100),
        
        sliderInput("seed_kernel_shape", 
                    "Seed kernel shape:", 
                    value = 1,
                    min = 0.1, 
                    max = 3,
                    step = 0.1)
      ), # End seed kernel panel
      
      wellPanel(
      
        h4("Pollen dispersal kernel options"),
        sliderInput("pollen_kernel_scale", 
                    "Pollen kernel scale:", 
                    value = 10,
                    min = 1, 
                    max = 100),
        sliderInput("pollen_kernel_shape", 
                    "Pollen kernel shape:", 
                    value = 1,
                    min = 0.1, 
                    max = 3,
                    step = 0.1)
      ),
      
      wellPanel(
        h4("Run simulation"),
        actionButton("initialize_new_simulation", "Initialize New Simulation"),
       sliderInput("n_steps", 
                   "Number of steps to run:", 
                   value = 1,
                   min = 1, 
                   max = 50),
       
       radioButtons("plotting_color", label = h3("Plotting color"),
                    choices = list("Maternal lineage" = "maternal", 
                                   "Paternal lineage" = "paternal"), 
                    selected = "maternal")
      ) # End simulation panel
         
  ), # End column of kernel options
  
  column(6, 
         h4("Plot of seed dispersal kernel"),
         plotOutput("seed_kernel_plot", height = "300px"),      
         
          h4("Plot of pollen dispersal kernel"),
                plotOutput("pollen_kernel_plot", height = "300px"),
         
         h4("Plot of simulation landscape"),
         plotOutput("sim_landscape_plot"),
         
         h4("Plot of simulation summary"),
         plotOutput("sim_summary_plot")
         )
  
  
  ) # End fluid row one
  
# fluidRow(
#   
#   column(3,
#          h4("Run simulation"),
#          sliderInput("n_steps", 
#                      "Number of steps to run:", 
#                      value = 1,
#                      min = 1, 
#                      max = 50),
#          
#          actionButton("run_new_simulation", "Run New Simulation"),
#          actionButton("update_plot", "Update plot of simulation")
# 
#   ),
#   column(9,
#          h4("Plot of simulation landscape"),
#          plotOutput("sim_summary_plot")
#   )
#   
#   
# ), # End row 4 
# 
# 
# fluidRow(
#   
#   column(12,
#          plotOutput("sim_landscape_plot")
#   )
#   
#   
# ) # End row 5
# 
# 

) # End fluid page
) # End Shiny ui