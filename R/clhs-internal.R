#' Objective function for cLHS
#'
#' @author Pierre Roudier
#' 
#' @keywords internal
#' 
#' @importFrom stats cor
#' @importFrom graphics hist
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

  if (n_cont_variables > 0) {
    
    cont_data_strata <- lapply(
      1:n_cont_variables, 
      function(i) list(data_continuous_sampled[, i, drop = TRUE], continuous_strata[, i, drop = TRUE])
    )
    
    cont_obj_sampled <- lapply(
      cont_data_strata, 
      function(x) hist(x[[1]], breaks = x[[2]], plot = FALSE)$counts
    )
    
    cont_obj_sampled <- matrix(unlist(cont_obj_sampled), ncol = n_cont_variables, byrow = FALSE)
  
    delta_obj_continuous <- rowSums(abs(cont_obj_sampled - eta))
    
    # Correlation of continuous data
    #
    cor_sampled <- suppressWarnings(cor(data_continuous_sampled))
    
    # when there's only one observation, cor() throws NAs - we set these to 1
    cor_sampled[is.na(cor_sampled)] <- 1 
    
    delta_obj_cor <- sum(abs(cor_mat - cor_sampled))
  
    } else {
    delta_obj_continuous <- 0
    delta_obj_cor <- 0
  }
  
  # Factor variables
  #
  n_factor_variables <- ncol(data_factor_sampled)
  
  if (n_factor_variables > 0) {
    factor_obj_sampled <- lapply(
      1:n_factor_variables, 
      function(x) {
        # Remove non-populated factor levels
        if (is.factor(data_factor_sampled[, x])) {
          data_factor_sampled[, x] <- droplevels(data_factor_sampled[, x])
        }
        
        res <- table(data_factor_sampled[, x])/nrow(data_factor_sampled)
        
        res
      }
    )
    
    delta_obj_factor <- lapply(
      1:n_factor_variables, 
      function(x) sum(abs(factor_obj_sampled[[x]] - factor_obj[[x]]))
    )

    # delta_obj_factor <- unlist(delta_obj_factor)#/length(delta_obj_factor) # do we need to ponder w/ the number of factors?
    delta_obj_factor <- abs( unlist(factor_obj_sampled) - unlist(factor_obj) )
    
  } else {
    delta_obj_factor <- 0
  }

  # Objective function
  #
  obj <- weights[[1]] * sum(delta_obj_continuous) + 
    weights[[2]] * sum(delta_obj_factor) + 
    weights[[3]] * delta_obj_cor

  # Returning results
  #
  list(
    obj = obj, 
    delta_obj_continuous = delta_obj_continuous, 
    delta_obj_factor = delta_obj_factor, 
    delta_obj_cor = delta_obj_cor
  )
}
