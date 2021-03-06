#' Objective function for cLHS
#'
#' @author Pierre Roudier
#' 
#' @keywords internal
#' 
#' @importFrom stats cor
#' @importFrom graphics hist
#' @noRd
#' 
.lhs_obj <- function(
  size,
  data_continuous_sampled,
  data_factor_sampled,
  continuous_strata,
  weights = list(numeric = 1, factor = 1, corelation = 1),
  cor_mat,
  factor_obj,
  eta = 1
  ) {

  # Continuous variables
  #
  n_cont_variables <- ncol(data_continuous_sampled)

  cont_data_strata <- lapply(
    1:n_cont_variables, 
    function(i) list(data = data_continuous_sampled[, i, drop = TRUE], strata = continuous_strata[, i, drop = TRUE]) 
  )
  
  cont_obj_sampled <- lapply(
    cont_data_strata, 
    function(x) hist(x[[1]], breaks = x[[2]], plot = FALSE)$counts
  )
  
  cont_obj_sampled <- matrix(unlist(cont_obj_sampled), ncol = n_cont_variables, byrow = FALSE)

  delta_obj_continuous <- rowSums(abs(cont_obj_sampled - eta))

  ## cont_obj_sampled tells use which histogram bucket is most out of whack with
  ## the target eta for each continuous covariate, adding across rows gets the
  ## per-bucket today for buckets 1 to n_samples across all covariates
  
  ## Compute information about how much each sample contributes to each histogram
  ## bucket across all features
  covariate.weights <- lapply(
    cont_data_strata, 
    function(x) {
      hst <- hist(
        x = x$data, 
        breaks = unique(x$strata), 
        plot = FALSE
      )
      res <- hst$counts[
        cut(x$data, breaks = unique(x$strata), labels = FALSE, include.lowest = TRUE)
      ]
      return(res)
    }
  )
  
  covariate.weights <- matrix(unlist(covariate.weights), ncol = n_cont_variables, byrow = FALSE)
  sample.weights <- rowSums(covariate.weights - eta)

  # Factor variables
  #
  n_factor_variables <- ncol(data_factor_sampled)
  if (n_factor_variables > 0) {
    factor_obj_sampled <- lapply(
      1:n_factor_variables, 
      function(x) {
        table(
          # Converting to factor to ensure all the levels in `factor obj` are
          # accounted for -- otherwise table drops the levels with "0" counts 
          factor(data_factor_sampled[, x], levels = names(factor_obj[[x]]))
        ) / nrow(data_factor_sampled)
    })
    
    delta_obj_factor <- lapply(
      1:n_factor_variables, 
      function(x) sum(abs(factor_obj_sampled[[x]] - factor_obj[[x]]))
    )

    delta_obj_factor <- unlist(delta_obj_factor)#/length(delta_obj_factor) # do we need to ponder w/ the number of factors?
    covariate.weights.factors <- lapply(
      1:n_factor_variables, 
      function(x) abs(factor_obj_sampled[[x]][data_factor_sampled[, x, drop = TRUE]] - factor_obj[[x]][data_factor_sampled[, x, drop = TRUE]])
    )
    
    covariate.weights.factors <- matrix(unlist(covariate.weights.factors), ncol = n_factor_variables, byrow = FALSE)
    
    sample.weights.factors <- rowSums(covariate.weights.factors)
    sample.weights.factors <- sample.weights.factors * nrow(data_factor_sampled) ## to put them on scale with continuous
    
  }
  else {
    delta_obj_factor <- 0
    sample.weights.factors <- rep(0, length(sample.weights))
  }
  
  # Correlation of continuous data
  #
  cor_sampled <- suppressWarnings(cor(data_continuous_sampled))
  
  # when there's only one observation, cor() throws NAs - we set these to 1
  cor_sampled[is.na(cor_sampled)] <- 1 

  delta_obj_cor <- sum(abs(cor_mat - cor_sampled))

  # Objective function
  #
  obj <- weights[[1]]*sum(delta_obj_continuous) + weights[[2]]*sum(delta_obj_factor) + weights[[3]]*delta_obj_cor

  # Returning results
  #
  list(
    obj = obj, 
    delta_obj_continuous = delta_obj_continuous, 
    delta_obj_factor = delta_obj_factor, 
    delta_obj_cor = delta_obj_cor, 
    # sample.weights = sample.weights
    sample.weights = weights[[1]] * sample.weights + weights[[2]] * sample.weights.factors
  )
}
