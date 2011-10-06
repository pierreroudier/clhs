## Conditioned latin Hypercube Sampling, after Minasny and McBratney, 2006.
##
## (c) Pierre Roudier, Landcare Research, 2011.
##

clhs <- function(
  x, # Continuous data
  size, # Number of samples you want
  n = 10000, # Number of max iterations
  tdecrease = 0.95,
  w1 = 1, # weight for continuous data
  w2 = 1, # weight for corelation among data
  w3 = 1, # weight for object data
  progress = TRUE # progress bar
  ) {

  # Detection of any factor data
  i_factor <- which(!sapply(x, is.numeric))
  n_factor <- length(i_factor)
  if (n_factor > 0) {
    data_continuous <- x[, -1*i_factor, drop = FALSE]
    data_factor <- x[, i_factor, drop = FALSE]
    # Creating a list storing the levels of each factor
    factor_levels <- apply(data_factor, 2, function(x) {
      ifelse(is.factor(x), res <- levels(x), res <-levels(factor(x)))
      res}
    )
  } else {
    data_continuous <- x
  }

  # annealing schedule
  temp <- 1
  metropolis <- 1

  n_data <- nrow(data_continuous) # Number of individuals in the data set

  # the edge of the strata
  xedge <- apply(data_continuous, 2, function(x) quantile(x, probs = seq(0, 1, length.out = size + 1)))

  # data correlation
  cor_mat <- cor(data_continuous)

  # for object/class variable
  if (n_factor == 0) {
    data_factor_sampled <- NULL
  }
  else {
    cobj <- apply(data_factor, 2, function(x) table(x)/n_data)
  }

  # initialise, pick randomly
  n_remainings <- n_data - size # number of individuals remaining unsampled
  i_sampled <- sample(1:n_data, size = size, replace = FALSE) # individuals randomly chosen
  i_unsampled <- setdiff(1:n_data, i_sampled) # individuals remaining unsampled
  data_continuous_sampled <- data_continuous[i_sampled, ] # sampled continuous data

  if (n_factor > 0)
    data_factor_sampled <- data_factor[i_sampled, ] # sampled factor data

  # objective function
  res <- .lhs_obj(size = size, data_continuous_sampled = data_continuous_sampled, data_factor_sampled = data_factor_sampled, xedge = xedge, cor_mat = cor_mat, cobj = cobj)

  obj <- res$obj # value of the objective function
#   isam <- res$isam
  dif <- res$dif #
#   iobj <- res$iobj

  # vector storing the values of the objective function
  obj_values <- vector(mode = 'numeric', length = n)

  # progress bar
  if (progress)
    pb <- txtProgressBar(min = 1, max = n, style = 3)

  for (i in 1:n) {

    # storing current values
    current <- list()
    current$obj <- obj
    current$i_sampled <- i_sampled
    current$i_unsampled <- i_unsampled
    current$dif <- dif

    if (runif(1) < 0.5) {
      # pick a random sample & swap with reservoir
      idx_unsampled <- sample(1:n_remainings, size = 1)
      idx_sampled <- sample(1:size, size = 1)
      # Retrieving indices values
      spl_unsampled <- i_unsampled[idx_unsampled]
      spl_sampled <- i_sampled[idx_sampled]
      # Swap these:
      i_sampled[idx_sampled] <- spl_unsampled
      i_unsampled[idx_unsampled] <- spl_sampled

      # creating new data sampled
      data_continuous_sampled[idx_sampled, ] <- data_continuous[idx_unsampled, ]
      if (n_factor > 0) {
        data_factor_sampled[idx_sampled, ] <- data_factor[idx_unsampled, ]
      }
    }
    else {
      # remove the worse sampled & resample
      worse <- max(dif)
      i_worse <- which(dif == worse)
      n_worse <- length(i_worse)

      # swap with reservoir
      spl_removed <- i_sampled[i_worse] # will be removed from the sampled set
      idx_added <- sample(1:n_remainings, size = n_worse) # will take their place
      i_sampled[i_worse] <- i_unsampled[idx_added] # replacing worst sampled by new pick
      i_unsampled[1:n_worse] <- spl_removed # replacing the worst pick in the reservoir

      # creating new data sampled
      data_continuous_sampled[i_worse, ] <- data_continuous[idx_added, ]
      if (n_factor > 0) {
        data_factor_sampled[i_worse, ] <- data_factor[idx_added, ]
      }
    }

    # calc obj
    res <- .lhs_obj(size = size, data_continuous_sampled = data_continuous_sampled, data_factor_sampled = data_factor_sampled, xedge = xedge, cor_mat = cor_mat, cobj = cobj)

    obj <- res$obj
#     isam <- res$isam
    dif <- res$dif
#     iobj <- res$iobj

    # Compare with previous iterations
    delta_obj <- obj - current$obj
    metropolis <- exp(-1*delta_obj/temp) + runif(1)*temp

    # If the optimum has been reached
    if (obj == 0) {
      obj_values[i] <- obj
      break
    }

    # Revert change
    if (delta_obj > 0 & runif(1) >= metropolis) {
#     if (delta_obj <=0 | runif(1) < metropolis) {# accept change
#      }
#     else {
      i_sampled <- current$i_sampled
      i_unsampled <- current$i_unsampled
      data_continuous_sampled <- data_continuous[i_sampled, ]
      if (n_factor > 0) {
        data_factor_sampled <- data_factor[i_sampled, ]
      }
      obj <- current$obj
      dif <- current$dif
    }

    # Storing the objective function value of the current iteration
    obj_values[i] <- obj

    # Temperature decrease
    if ((i %% 10) == 0) {
      temp <- temp*tdecrease
    }

    # Update progress bar
    if (progress)
      setTxtProgressBar(pb, i)
  }

  # Close progress bar
  if (progress)
    close(pb)

  res <- list(index_samples = i_sampled, sampled_data = data_continuous_sampled, obj_function = obj_values)
  class(res) <- "cLHS_result"
  res
}

## objective function for LHS

.lhs_obj <- function(
  size,
  data_continuous_sampled,
  data_factor_sampled,
  xedge,
  w1 = 1,
  w2 = 1,
  cor_mat,
  w3 = 1,#,
#   factor_levels,
  cobj
  ) {

  # Continuous variables
  n_cont_variables <- ncol(data_continuous_sampled)

  cdat <- lapply(1:n_cont_variables, function(i) list(data_continuous_sampled[, i], xedge[, i]) )
  isam <- lapply(cdat, function(x) hist(x[[1]], breaks = x[[2]], plot = FALSE)$counts)
  isam <- matrix(unlist(isam), ncol = n_cont_variables, byrow = FALSE)

  dif <- rowSums(abs(isam - 1))

  # Correlation of continuous data
  cor_sampled <- cor(data_continuous_sampled)
  dc <- sum(sum(abs(cor_mat - cor_sampled)))

  # Factor variables
  n_factor_variables <- ncol(data_factor_sampled)
  iobj <- lapply(1:n_factor_variables, function(x) table(data_factor_sampled[,x])/nrow(data_factor_sampled) - cobj[[x]])
  do <- lapply(1:n_factor_variables, function(x) sum(abs(iobj[[x]] - cobj[[x]])))
  do <- sum(unlist(do)/length(do))

  obj <- w1*sum(dif) + w2*dc + w3*do

  list(obj = obj, isam = isam, dif = dif)
}

## Basic plotting function

plot.cLHS_result <- function(x, ...) {
  plot(1:length(x$obj_function), x$obj_function, type = 'l', xlab = "Iteration", ylab = "Energy function", ...)
}