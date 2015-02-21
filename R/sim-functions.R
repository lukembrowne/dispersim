

## Move forward x generations
stepGeneration <- function(world, steps = 1){
  while(steps > 0){
    
    world$step_counter <- world$step_counter + 1
    world <- survival(world)
    world <- reproduce(world)
    world <- saveIdAlive(world)
    #generate summary
    # print status
    steps <- steps - 1 

      
  } 
  return(world)
}        

saveIdAlive <- function(world){
  world$registry[[world$step_counter]] <- list(
    adults_alive = which(world$data$type == "adult" & world$data$alive == TRUE),
    seedlings_alive = which(world$data$type == "seedling" & world$data$alive == TRUE))
  return(world)
}