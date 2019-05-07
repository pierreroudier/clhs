#' @include clhs-internal.R
#' @rdname clhs
#' @method clhs data.frame
#' @importFrom stats quantile runif 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
clhs.data.frame <- function(
  x, # data.frame
  size, # Number of samples you want
  include = NULL, # row index of data that must be in the final sample
  cost = NULL, # Number or name of the attribute used as a cost
  iter = 10000, # Number of max iterations
  temp = 1, # initial temperature
  tdecrease = 0.95, # temperature decrease rate
  weights = list(numeric = 1, factor = 1, correlation = 1), # weight for continuous data , weight for correlation among data, weight for object data
  eta = 1,
  obj.limit = -Inf, # Stopping criterion
  length.cycle = 10, # Number of cycles done at each constant temperature value
  simple = TRUE, # only return selected indices (if false, return a more complex S3 object)
  progress = TRUE, # progress bar,
  track = NULL # just to have the cost computed without having it guiding the process
) {
  
  # Temperature decrease rate should be < 1
  if (tdecrease >= 1) stop("tdecrease should be < 1")
  
  ## No cost taking into account during optimisation
  if (is.null(cost)) {
    cost_mode <- FALSE
    
    if (!is.null(track)) {
      ## Track mode: tracking cost (without taking it into account during optimisation)
      
      # Get the id of the column used as a cost attribute
      if (is.numeric(track)) i_cost <- track
      else i_cost <- which(names(x) == track)
      
      # Get cost attribute column 
      cost <- x[ , i_cost, drop = FALSE]
      # Remove cost attribute from attribute table
      x <- x[, -1*i_cost, drop = FALSE]
      
      # Flags
      cost_mode <- TRUE
      track_mode <- TRUE
      
    } else { 
      ## No cost tracking, just plain optimisation of attributes
      track_mode <- FALSE
    }
    
  } else {
    ## Cost is taken into account for optimisation
    
    # Get the id of the column used as a cost attribute
    if (is.numeric(cost)) i_cost <- cost
    else i_cost <- which(names(x) == cost)
    
    if (!length(i_cost)) stop("Could not find the cost attribute.") 
    
    # Get cost attribute column 
    cost <- x[ , i_cost, drop = FALSE]
    # Remove cost attribute from attribute table
    x <- x[, -1*i_cost, drop = FALSE]
    
    # Si include, cost is 0
    if (!is.null(include)) {
      cost[include, ] <- 0
    }
    
    # Flags
    cost_mode <- TRUE
    track_mode <- FALSE # cost is taken into account, therefore computed
    
  }
  
  if (!is.null(include)) {
    if (size <= length(include)) {
      stop(paste0("size (", size, ") should be larger than length of include (", length(include), ")"))
    }
  }
  
  # Detection of any factor data
  i_factor <- which(!sapply(x, is.numeric))
  n_factor <- length(i_factor)
  n_continuous <- ncol(x) - n_factor
    
  # If there is a mix of continuous and categorical data
  if (n_factor > 0 & n_continuous > 0) { 
    data_continuous <- x[, -1*i_factor, drop = FALSE]
    data_factor <- x[, i_factor, drop = FALSE]
    # Creating a list storing the levels of each factor
    factor_levels <- apply(data_factor, 2, function(x) {
      ifelse(is.factor(x), res <- levels(x), res <- levels(factor(x)))
      res}
    )
  } else if (n_factor == 0 & n_continuous > 0) {
    # Continuous data but no categorical data
    data_continuous <- x
  } else if (n_factor > 0 & n_continuous == 0) {
    # categorical data but no continuous data
    data_factor <- x
  } else {
    stop("Wrong types in input data.frame.")
  }
  
  metropolis <- exp(-1*0/temp) # Initial Metropolis value
  n_data <- nrow(x) # Number of individuals in the data set
  
  if (n_continuous > 0) {
    # Edge of the strata
    continuous_strata <- apply(
      data_continuous, 
      2, 
      function(x) {
        quantile(x, probs = seq(0, 1, length.out = size + 1), na.rm = TRUE)
      }
    )
    
    # Data correlation
    cor_mat <- cor(data_continuous, use = "complete.obs")
  } else {
    continuous_strata <- cor_mat <- NULL
  }
  
  # For object/class variable
  if (n_factor == 0) data_factor_sampled <- data.frame()
  else factor_obj <- alply(data_factor, 2, function(x) table(x)/n_data)
  
  # Mandatory data in the sample
  sampled_size <- size - length(include)
  not_included <- setdiff(1:n_data, include)
  
  # initialise, pick randomly
  n_remainings <- n_data - size # number of individuals remaining unsampled
  i_sampled <- c(sample(not_included, size = sampled_size, replace = FALSE), include) # individuals randomly chosen
  i_unsampled <- setdiff(1:n_data, i_sampled) # individuals remaining unsampled
  
  if (n_continuous > 0) {
    data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE] # sampled continuous data
  } else {
    data_continuous_sampled <- data.frame()
  }
  if (n_factor > 0) {
    data_factor_sampled <- data_factor[i_sampled, , drop = FALSE] # sampled factor data
  } else {
    data_factor_sampled <- data.frame()
  }
  
  # objective function
  res <- .lhs_obj(
    size = size, 
    data_continuous_sampled = data_continuous_sampled, 
    data_factor_sampled = data_factor_sampled, 
    continuous_strata = continuous_strata, 
    cor_mat = cor_mat, 
    factor_obj = factor_obj, 
    weights = weights, 
    eta = eta
  )
  
  obj <- res$obj # value of the objective function
  delta_obj_continuous <- res$delta_obj_continuous
  delta_obj_factor <- res$delta_obj_factor
  
  if (cost_mode) {
    # (initial) operational cost
    op_cost <- sum(cost[i_sampled, ])
    # vector storing operational costs
    op_cost_values <- vector(mode = 'numeric', length = iter)
  } else op_cost_values <- NULL
  
  # vector storing the values of the objective function
  obj_values <- vector(mode = 'numeric', length = iter)
  
  # progress bar
  if (progress) pb <- txtProgressBar(min = 1, max = iter, style = 3)
  flip <- -1
  
  for (i in 1:iter) {
    print(paste0(i, " | ", flip < 0.5, "|", length(i_sampled)))
    
    # storing previous values
    previous <- list()
    previous$obj <- obj
    previous$i_sampled <- i_sampled
    previous$i_unsampled <- i_unsampled
    previous$delta_obj_continuous <- delta_obj_continuous
    previous$delta_obj_factor <- delta_obj_factor
    
    if (cost_mode) previous$op_cost <- op_cost
    
    flip <- runif(1)
    
    if (flip < 0.5) {
      # pick a random sampled point and random unsampled point
      # 
      rm_candidates <- setdiff(i_sampled, include)
      spl_removed <- sample(rm_candidates, size = 1, replace = FALSE) 
      spl_added <- sample(i_unsampled, size = 1, replace = FALSE)
      
      # swap the added and removed points from i_sampled and i_unsampled
      #
      
      # 1. remove sample from i_sampled 
      i_sampled <- setdiff(i_sampled, c(include, spl_removed))
      
      # 2. replace it by new sampled point 
      i_sampled <- c(i_sampled, spl_added, include)
      
      # 3. remove the added sampled point from the pool of unsampled points
      i_unsampled <- setdiff(i_unsampled, spl_added)
      
      # 4. add the removed sampled back into pool of unsampled points
      i_unsampled <- c(i_unsampled, spl_removed)
      
      if (any(duplicated(i_sampled))) browser() 
      if (any(i_sampled %in% i_unsampled)) browser()
      
      # creating new data sampled
      if (n_continuous > 0) data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE]
      if (n_factor > 0) data_factor_sampled <- data_factor[i_sampled, , drop = FALSE]
      
    }
    else {
      
      browser()
      
      # Compute objective function for each candidate
      individual_obj <- weights[[1]] * res$delta_obj_continuous + 
        weights[[2]] * res$delta_obj_factor + 
        weights[[3]] * res$delta_obj_cor
      
      # remove the worse sampled & resample
      
      # worse <- max(delta_obj_continuous[!i_sampled %in% include])
      # i_worse <- which(delta_obj_continuous[!i_sampled %in% include] == worse)
      # 
      worse <- max(individual_obj[!i_sampled %in% include])
      i_worse <- which(individual_obj[!i_sampled %in% include] == worse)
      
      # If there's more than one worse candidate, we pick one at random
      if (length(i_worse) > 1) i_worse <- sample(i_worse, size = 1)
      
      # swap with reservoir
      
      spl_removed <- setdiff(i_sampled, include)[i_worse]
      spl_added <- sample(i_unsampled, size = 1, replace = FALSE)
      
      # spl_removed <- setdiff(i_sampled, include)[i_worse] # will be removed from the sampled set.
      # idx_added <- sample(1:n_remainings, size = 1, replace = FALSE) # new candidate that will take their place
      # i_sampled <- setdiff(i_sampled, include)[-i_worse]
      
      # if (any(is.na(i_sampled))) browser()
      # 
      # i_sampled <- c(i_sampled, i_unsampled[idx_added], include)
      # i_unsampled <- i_unsampled[-idx_added]
      # i_unsampled <- c(i_unsampled, spl_removed)
      
      # swap the added and removed points from i_sampled and i_unsampled
      #
      
      # 1. remove sample from i_sampled 
      i_sampled <- setdiff(i_sampled, c(include, spl_removed))
      
      # 2. replace it by new sampled point 
      i_sampled <- c(i_sampled, spl_added, include)
      
      # 3. remove the added sampled point from the pool of unsampled points
      i_unsampled <- setdiff(i_unsampled, spl_added)
      
      # 4. add the removed sampled back into pool of unsampled points
      i_unsampled <- c(i_unsampled, spl_removed)
      
      if (any(i_sampled %in% i_unsampled)) browser()
      
      if (any(duplicated(i_sampled))) stop("duplicated point(s) in samples")
      if (any(i_sampled %in% i_unsampled)) stop("point(s) in both sampled and unsampled datasets")
      
      # creating new data sampled
      if (n_continuous > 0) data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE]
      if (n_factor > 0) data_factor_sampled <- data_factor[i_sampled, , drop = FALSE]
    }
    
    # calc obj
    res <- .lhs_obj(
      size = size, 
      data_continuous_sampled = data_continuous_sampled, 
      data_factor_sampled = data_factor_sampled, 
      continuous_strata = continuous_strata, 
      cor_mat = cor_mat, 
      factor_obj = factor_obj, 
      weights = weights, 
      eta = eta
    )
    
    obj <- res$obj
    delta_obj_continuous <- res$delta_obj_continuous
    delta_obj_factor <- res$delta_obj_factor
    
    # Compare with previous iterations
    delta_obj <- obj - previous$obj
    metropolis <- exp(-1*delta_obj/temp) #+ runif(1)*temp
    
    if (cost_mode) {
      # op costs
      op_cost <- sum(cost[i_sampled, ])
      delta_cost <- op_cost - previous$op_cost
      if (track_mode) metropolis_cost <- Inf # runif(1) >= Inf is always FALSE
      else metropolis_cost <- exp(-1*delta_cost/temp)
    }
    else metropolis_cost <- Inf # runif(1) >= Inf is always FALSE
    
    # If the optimum has been reached
    if (obj <= obj.limit) {
      warning("\nThe objective function has reached its minimum value, as specified by the obj.limit option.")
      
      if (progress) {
        setTxtProgressBar(pb, i)
        close(pb)
      }
      
      obj_values[i] <- obj
      
      if (cost_mode) op_cost_values[i] <- op_cost
      
      break
    }
    
    # Revert change
    if (delta_obj > 0 & runif(1) >= metropolis | runif(1) >= metropolis_cost) {
      i_sampled <- previous$i_sampled
      i_unsampled <- previous$i_unsampled
      
      if (n_continuous > 0) data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE]
      if (n_factor > 0) data_factor_sampled <- data_factor[i_sampled, , drop = FALSE]
      
      obj <- previous$obj
      delta_obj_continuous <- previous$delta_obj_continuous
      delta_obj_factor <- previous$delta_obj_factor
      
      if (cost_mode) op_cost <- previous$op_cost
    }
    
    # Storing the objective function value of the current iteration
    obj_values[i] <- obj
    
    if (cost_mode) op_cost_values[i] <- op_cost
    
    # Temperature decrease
    if ((i %% length.cycle) == 0) temp <- temp*tdecrease
    
    # Update progress bar
    if (progress) setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  if (progress) close(pb)
  
  if (n_factor > 0 & n_continuous > 0) {
    sampled_data <- data.frame(data_continuous_sampled, data_factor_sampled)
    # reordering cols
    sampled_data <- sampled_data[, names(x)]
  } else if(n_factor == 0 & n_continuous > 0) {
    sampled_data <- data_continuous_sampled 
  } else if(n_factor > 0 & n_continuous == 0) {
    sampled_data <- data_factor_sampled 
  } else {
    stop("Wrong types in input data.frame.")
  }
  
  # Simple output - just the sampled object
  if (simple) res <- i_sampled
  else {
    # Making up the object to be returned
    res <- list(
      initial_object = x,
      index_samples = i_sampled, 
      sampled_data = sampled_data, 
      obj = obj_values,
      cost = op_cost_values
    )
    class(res) = c("cLHS_result","list")
  }
  
  res
}