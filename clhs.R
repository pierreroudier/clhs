## Conditioned latin Hypercube Sampling, after Minasny and McBratney, 2006.
##
## (c) Pierre Roudier, Landcare Research, 2011.
##

clhs <- function(
  x, # Continuous data
  size, # Number of samples you want
  niter = 10000, # Number of max iterations
  w1 = 1, # weight for continuous data
  w2 = 1, # weight for corelation among data
  w3 = 1 # weight for object data
) {

  # Factor data
  i_factor <- which(sapply(x, class) %in% c("factor", "character"))
  n_factor <- length(i_factor)
  if (n_factor > 0) {
    xobj <- x[, i_factor, drop = FALSE]
    vobj <- levels(xobj)
  }

  # annealing schedule
  temp <- 1
  tdecrease <- 0.95
  metrop <- 1

  n_data <- nrow(x)
  n_variables <- ncol(x)

  # the edge of the strata
#   step <- 100/size
#   Pl <- seq(0, 100, by = step)
#   pc <- (seq(0.5, n_data) - 0.5) / n_data

  xedge <- matrix(NA, ncol = n_variables, nrow = size + 1)
  for (j in 1:n_variables){
    xedge[, j] <- quantile(x[, j], probs = seq(0, 1, length.out = size + 1))
#     Y <- sort(x[,j])
#     I <- order(x[,j])
#       xquant(I,j)=pc';
  }

#   browser()

  # target value
#   ibest <- matrix(1, ncol = n_variables, nrow = size)

  # data correlation
  corr <- cor(x)

  # for object/class variable
  if (n_factor == 0) {
      nobj <- 0
      cobj <- 0
      w3 <- 0
      xobs <- NULL
  }
  else {
      nobj <- length(vobj)
      # For each class value (?)
#       for (j in 1:nobj) {
#         cobj[j] <- sum(xobj == vobj[j])
#       }
      cobj <- table(xobj)
      cobj <- (cobj/n_data)*size # target value
  }

  # initialise, pick randomly
  nunsam <- n_data - size
  idat <- randperm(1:n_data)
  ipick <- idat[1 : size]
  iunsam <- idat[(size + 1) : n_data]
  xsam <- x[ipick, ]
  if (nobj > 0)
    xobs <- xobj[ipick, ]


  # objective function
  res <- lhs_obj(size = size, n_variables = n_variables, xsam = xsam, xobs = xobs, xedge = xedge, corr = corr)
  obj <- res$obj
  isam <- res$isam
  dif <- res$dif
#   iobj <- res$iobj

  obj_values <- list()

  for (il in 1:niter) {
    idx <- randperm(1:nunsam)
    iunsam <- iunsam[idx]

    # storing stuff
    objsave <- obj
    isave <- ipick
    iusave <- iunsam
    dsave <- dif

    if (runif(1) < 0.5) {
      # pick a random sample & swap with reservoir
      iw <- randperm(1:size)
      ich <- iunsam[1]
      jch <- ipick[iw[1]]
      ipick[iw[1]] <- ich
      iunsam[1] <- jch
      xsam[iw[1], ] <- x[ich, ]
      if (nobj > 0) {
        xobs[iw[1], ] <- xobj[ich, ]
      }
    }
    else {
      # remove the worse sampled & resample
      worse <- max(dif)
      iworse <- which(dif == worse)
      nworse <- length(iworse)

      # swap with reservoir
      jch <- ipick[iworse]
      kch <- iunsam[1:nworse]
      ipick[iworse] <- kch
      iunsam[1:nworse] <- jch
      xsam[iworse, ] <- x[kch, ]
      if (nobj > 0) {
        xobs[iworse, ] <- xobj[kch, ]
      }
    }

    # calc obj
    res <- lhs_obj(size = size, n_variables = n_variables, xsam = xsam, xobs = xobs, xedge = xedge, corr = corr)
    obj <- res$obj
    isam <- res$isam
    dif <- res$dif
#     iobj <- res$iobj

    #  compare with previous iterations
    de <- obj - objsave
    metrop <- exp(-1*de/temp) + runif(1)*temp

    if (obj == 0) {
      obj_values[[il]] <- obj
      break
    }

    if (de <=0 | runif(1) < metrop) {# accept change
    }
    else {
      # revert, no changes
      ipick <- isave
      iunsam <- iusave
      xsam <- x[ipick, ]
      if (nobj > 0) {
        xobs <- xobj[ipick, ]
      }
      obj <- objsave
      dif <- dsave
    }

    obj_values[[il]] <- obj
    if ((il %% 10) ==0) {
      temp <- temp*tdecrease
    }

    # calc the final obj function
    res <- lhs_obj(size = size, n_variables = n_variables, xsam = xsam, xobs = xobs, xedge = xedge, corr = corr)
    obj <- res$obj
    isam <- res$isam
    dif <- res$dif
#     iobj <- res$iobj
  }
  list(ipick, xsam, obj_values, isam)
}

## Random permutations
randperm <- function(x, ...) {
  sample(x, size = length(x), replace = FALSE, ...)
}

## objective function for LHS

lhs_obj <- function(
  size,
  n_variables,
  xsam,
  xobs,
  xedge,
  w1 = 1,
  w2 = 1,
  corr,
  w3 = 1#,
#   vobj,
#   cobj
  ) {

  # data in quantiles
  hsam <- isam <- matrix(NA, ncol = n_variables, nrow = size)
  for (j in 1:n_variables) {
    hsam[, j] <- hist(xsam[, j], breaks = xedge[, j], plot = FALSE)$counts
    isam[, j] <- hsam[1:size, j]
  }
  dif <- sum(abs(isam - 1))

  # correlation
  csam <- cor(xsam)
  dc <- sum(sum(abs(corr - csam)))

#   # object data
#   nobj <- length(vobj)
#   for (j in 1:nobj) {
#     iobj[j] <- sum(xobs == vobj[j])
#   }
#   do <- sum(abs(iobj - cobj))

  obj <- w1*sum(dif) + w2*dc #+ w3*do

  list(obj = obj, isam = isam, dif = dif)
}
