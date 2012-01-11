# Produces a plot illustrating teh results of the cLHS sampling.
# - a plot of the objective function
# - histograms or density plots of the sampled attributes
#
plot.cLHS_result <- function(
  x,
  modes = "obj",
  ...
  ){

  require(ggplot2)

  # Number of canvas to init
  n_views <- length(modes)
  iter <- 1:length(x$obj_function)
  obj <- x$obj_function
  df_obj <- data.frame(iter, obj)

  pl <- list()

  if ("obj" %in% modes) {
  objective_plot <- ggplot(df_obj) + geom_line(aes(x = iter, y = obj)) + labs(x = "Iteration", y = "Objective function") + theme_bw()  + opts(title = "Evolution of the objective function")

  pl[[length(pl) + 1]] <- objective_plot
  }

  if (any(c("hist", "dens" )%in% modes)) {
    
    if (all(c("hist", "dens") %in% modes))
      stop('"hist" and "dens" modes are mutually exclusive.')

    init <- x$initial_object
    spl <- x$sampled_data
    
    if (.is.spatial(init)) {
      if (inherits(init, "Spatial")) {
        init <- init@data
        spl <- spl@data
      }
      if (inherits(init, "Raster")) {
        init <- as.data.frame(init)
        spl <- as.data.frame(spl)
      }
    }

    # Separate continuous from factor variables
    i_factor <- which(!sapply(init, is.numeric))
    i_continuous <- setdiff(1:ncol(init), i_factor)
    n_factor <- length(i_factor)
    n_continuous <- length(i_continuous)

    if (n_factor > 0) {

    }

    init_continuous <- init[, i_continuous, drop = FALSE]
    spl_continuous <- spl[, i_continuous, drop = FALSE]
    init_factor <- init[, i_factor, drop = FALSE]
    spl_factor <- spl[, i_factor, drop = FALSE]

    # initiate an "id" column
    idcolname <- .create_unique_colname(id = "id", nm = names(init))
    
    if (n_factor > 0) {
      init_factor[[idcolname]] <- "init"
      spl_factor[[idcolname]] <- "spl"
      
      # merge df
      df_hist_factor <- melt(rbind(init_factor, spl_factor), idcolname)
      
      # Plot for factors (bar counts)
      distrib_factor <- ggplot(df_hist_factor) + geom_bar(aes_string(x = 'value', fill = idcolname), position = 'dodge') + facet_wrap(~variable, scales = "free") + theme_bw() + opts(title = "Discrete variables")

      pl[[length(pl) + 1]] <- distrib_factor
    }
    
    if (n_continuous > 0) {
      init_continuous[[idcolname]] <- "init"
      spl_continuous[[idcolname]] <- "spl"
  
      # merge df
      df_hist_continuous <- melt(rbind(init_continuous, spl_continuous), idcolname)
      
      # Plot for continuous
      if ('dens' %in% modes)
        distrib_continuous <- ggplot(df_hist_continuous) + geom_density(aes_string(x = 'value', fill = idcolname), alpha = 0.5) + facet_wrap(~ variable, scales = "free") + theme_bw()  + opts(title = "Continuous variables")
      if ('hist' %in% modes)
        distrib_continuous <- ggplot(df_hist_continuous) + geom_histogram(aes_string(x = 'value', fill = idcolname), position = 'dodge') + facet_wrap(~ variable, scales = "free") + theme_bw()  + opts(title = "Continuous variables")

      pl[[length(pl) + 1]] <- distrib_continuous
    }
  }

  if (length(pl) > 1) {
    grid.newpage()

    if (length(pl) == 2) {
      pushViewport(viewport(layout = grid.layout(1, 2)))
      vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

      print(pl[[1]], vp = vplayout(1, 1))
      print(pl[[2]], vp = vplayout(1, 2)) 
    }
    else {
      pushViewport(viewport(layout = grid.layout(2, 2)))
      vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

      print(pl[[1]], vp = vplayout(1:2, 1))
      print(pl[[2]], vp = vplayout(1, 2)) 
      print(pl[[3]], vp = vplayout(2, 2)) 
    }
  }
  else {
    print(pl[[1]])
  }

}