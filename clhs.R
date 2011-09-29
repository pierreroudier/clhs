## Conditioned latin Hypercube Sampling, after Minasny and McBratney, 2006.
##
## (c) Pierre Roudier, Landcare Research, 2011.
##

clhs <- function(
  x, # Continuous data
  size, # Number of samples you want
  n = 10000, # Number of max iterations
  w1 = 1, # weight for continuous data
  w2 = 1, # weight for corelation among data
  w3 = 1 # weight for object data
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
  tdecrease <- 0.95
  metrop <- 1

  n_data <- nrow(data_continuous) # Number of individuals in the data set
  n_cont_variables <- ncol(data_continuous) # Number of continuous variables

  # the edge of the strata
  xedge <- apply(data_continuous, 2, function(x) quantile(x, probs = seq(0, 1, length.out = size + 1)))

  # data correlation
  cor_mat <- cor(data_continuous)

  # for object/class variable
  if (n_factor == 0) {
    data_factor_sampled <- NULL
  }
  else {
      nobj <- length(factor_levels)
      # For each class value (?)
#       for (j in 1:nobj) {
#         cobj[j] <- sum(data_factor == factor_levels[j])
#       }
      cobj <- table(data_factor)
      cobj <- (cobj/n_data)*size # target value
  }

  # initialise, pick randomly
  n_remainings <- n_data - size # number of individuals remaining unsampled
  i_sampled <- sample(1:n_data, size = size, replace = FALSE) # individuals randomly chosen
  i_unsampled <- setdiff(1:n_data, i_sampled) # individuals remaining unsampled
  data_continuous_sampled <- data_continuous[i_sampled, ] # sampled continuous data

  if (n_factor > 0)
    data_factor_sampled <- data_factor[i_sampled, ] # sampled factor data

  # objective function
  res <- lhs_obj(size = size, n_cont_variables = n_cont_variables, data_continuous_sampled = data_continuous_sampled, data_factor_sampled = data_factor_sampled, xedge = xedge, cor_mat = cor_mat)

  obj <- res$obj # value of the objective function
#   isam <- res$isam
  dif <- res$dif #
#   iobj <- res$iobj

  # vector storing the values of the objective function
  obj_values <- vector(mode = 'numeric', length = n)

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
      # Retrieing indices values
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
      iworse <- which(dif == worse)
      nworse <- length(iworse)

      # swap with reservoir
      jch <- i_sampled[iworse] # will be removed
      kch <- i_unsampled[1:nworse] # will take their place
      i_sampled[iworse] <- kch # replacing worst sampled by new pick
      i_unsampled[1:nworse] <- jch # replacing the worst pick in the reservoir

      # creating new data sampled
      data_continuous_sampled[iworse, ] <- data_continuous[kch, ]
      if (n_factor > 0) {
        data_factor_sampled[iworse, ] <- data_factor[kch, ]
      }
    }

    # calc obj
    res <- lhs_obj(size = size, n_cont_variables = n_cont_variables, data_continuous_sampled = data_continuous_sampled, data_factor_sampled = data_factor_sampled, xedge = xedge, cor_mat = cor_mat)

    obj <- res$obj
#     isam <- res$isam
    dif <- res$dif
#     iobj <- res$iobj

    #  compare with previous iterations
    de <- obj - current$obj
    metrop <- exp(-1*de/temp) + runif(1)*temp

    if (obj == 0) {
      obj_values[i] <- obj
      break
    }

    if (de <=0 | runif(1) < metrop) {# accept change
    }
    else {
      # revert, no changes
      i_sampled <- current$i_sampled
      i_unsampled <- current$i_unsampled
      data_continuous_sampled <- data_continuous[i_sampled, ]
      if (n_factor > 0) {
        data_factor_sampled <- data_factor[i_sampled, ]
      }
      obj <- current$obj
      dif <- current$dif
    }

    obj_values[i] <- obj
    if ((i %% 10) ==0) {
      temp <- temp*tdecrease
    }

    # calc the final obj function
    res <- lhs_obj(size = size, n_cont_variables = n_cont_variables, data_continuous_sampled = data_continuous_sampled, data_factor_sampled = data_factor_sampled, xedge = xedge, cor_mat = cor_mat)
    obj <- res$obj
#     isam <- res$isam
    dif <- res$dif
#     iobj <- res$iobj
  }
  list(i_sampled, data_continuous_sampled, obj_values)
}

## Random permutations
randperm <- function(x, ...) {
  sample(x, size = length(x), replace = FALSE, ...)
}

## objective function for LHS

lhs_obj <- function(
  size,
  n_cont_variables,
  data_continuous_sampled,
  data_factor_sampled,
  xedge,
  w1 = 1,
  w2 = 1,
  cor_mat,
  w3 = 1#,
#   factor_levels,
#   cobj
  ) {

  # data in quanties
  hsam <- isam <- matrix(NA, ncol = n_cont_variables, nrow = size)
  for (j in 1:n_cont_variables) {
    hsam[, j] <- hist(data_continuous_sampled[, j], breaks = xedge[, j], plot = FALSE)$counts
    isam[, j] <- hsam[1:size, j]
  }

  dif <- sum(abs(isam - 1))

  # correlation of continuous data
  csam <- cor(data_continuous_sampled)
  dc <- sum(sum(abs(cor_mat - csam)))

#   # object data
#   nobj <- length(factor_levels)
#   for (j in 1:nobj) {
#     iobj[j] <- sum(data_factor_sampled == factor_levels[j])
#   }
#   do <- sum(abs(iobj - cobj))

  obj <- w1*dif + w2*dc #+ w3*do
  list(obj = obj, isam = isam, dif = dif)
}
