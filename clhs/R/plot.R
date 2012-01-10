# Produces a plot illustrating teh results of the cLHS sampling.
# - a plot of the objective function
# - histograms or density plots of the sampled attributes
#
plot.cLHS_result <- function(
  x,
  modes = "obj"
  ){

  require(lattice)

  # Number of canvas to init
  n_views <- length(modes)
  
  df <- data.frame(iteration = 1:length(x$obj_function), objective = x$obj_function)
  xyplot(objective ~ iteration, data = df, xlab = "Iteration", ylab = "Objective function", type = 'l') 
}