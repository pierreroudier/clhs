#' Conditioned Latin Hypercube Sampling
#' 
#' Implementation of the conditioned Latin hypercube sampling, as published by
#' Minasny and McBratney (2006) and the DLHS variant method (Minasny and 
#' McBratney, 2010). These methods propose to stratify sampling in
#' presence of ancillary data. An extension of this method, which propose to
#' associate a cost to each individual and take it into account during the
#' optimisation process, is also proposed (Roudier et al., 2012).
#' 
#' For the DLHS method, the original paper defines parameter \code{b} as the importance 
#' of the edge of the distributions. A matrix \code{eta} (size N x K, where N is the size of 
#' the final sample and K the number of continuous variables) is defined, to 
#' compute the objective function of the algorithm, where each column equal the 
#' vector (b, 1, ..., 1, b) in order to give the edge of the distribution a 
#' probability b times higher to be sampled. In our function, instead of define 
#' the \code{b} parameter, users can defined their own \code{eta} matrix so that they 
#' can give more complex probabilty design of sampling each strata of the 
#' distribution instead of just be able to give more importance to both edges of 
#' the distribution.
#' 
#' @name clhs
#' @aliases clhs clhs.data.frame clhs.SpatialPointsDataFrame clhs.Raster
#' clhs
#' 
#' @docType methods
#' 
#' @param x A \code{data.frame}, \code{SpatialPointsDataFrame} or \code{Raster}
#' object.
#' @param size A non-negative integer giving the total number of items to select
#' @param include A numeric vector giving the indices of the rows from \code{x} that must be 
#' included in the selected items. For the cost-constrained cLHS method, cost of 
#' these mandatory samples is set to 0. If NULL (default), all data are randomly 
#' choosen according to the classic cLHS method. If \code{include} is not NULL,
#' argument \code{size} must inlcude the total size of the final sample i.e. the
#' size of mandatory samples given by \code{include} plus the size of the randomly
#' chosen samples to pick.
#' @param cost A character giving the name or an integer giving the index of
#' the attribute in \code{x} that gives a cost that can be use to constrain the
#' cLHS sampling. If NULL (default), the cost-constrained implementation is not
#' used.
#' @param track A character giving the name or an integer giving the index
#' of the attribute in \code{x} that gives a cost associated with each
#' individual. However, this method will only track the cost - the sampling
#' prrocess will not be constrained by this attribute. If NULL (default), this
#' option is not used.
#' @param iter A positive number, giving the number of iterations for the
#' Metropolis-Hastings annealing process. Defaults to 10000.
#' @param temp The initial temperature at which the simulated annealing
#' begins. Defaults to 1.
#' @param tdecrease A number betwen 0 and 1, giving the rate at which
#' temperature decreases in the simulated annealing process. Defaults to 0.95.
#' @param weights A list a length 3, giving the relative weights for
#' continuous data, categorical data, and correlation between variables.
#' Defaults to \code{list(numeric = 1, factor = 1, correlation = 1)}.
#' @param eta Either a number equal 1 to perform a classic cLHS or a constrained 
#' cLHS or a matrix to perform a cLHS that samples more on the edge of the
#' distibutions (DLHS, see details)
#' @param obj.limit The minimal value at which the optimisation is stopped.
#' Defaults to \code{-Inf}.
#' @param length.cycle The duration (number of iterations) of the
#' isotemperature steps. Defaults to 10.
#' @param progress TRUE or FALSE, displays a progress bar.
#' @param simple TRUE or FALSE. If set to TRUE, only the indices of the
#' selected samples are returned, as a numeric vector. If set to FALSE, a
#' cLHS_result object is returned (takes more memory but allows to make use of
#' cLHS_results methods such as \code{plot.cLHS_result}).
#' @param ... additional parameters passed to \code{clhs}
#' 
#' @return * If the \code{simple} option is set to TRUE (default behaviour): A
#' numeric vector containing the indices of the selected samples is returned
#' 
#' * If the \code{simple} option is set to FALSE: An object of class
#' \code{cLHS_result}, with the following elements: \item{index_samples}{a
#' vector giving the indices of the chosen samples.} \item{sampled_data}{the
#' sampled data.frame.} \item{obj}{a vector giving the evolution of the
#' objective function throughout the Meropolis-Hastings iterations.}
#' \item{cost}{a vector giving the evolution of the cost function throughout
#' the Metropolis-Hastings iterations (if available).}
#' @author Pierre Roudier
#' @seealso \code{\link{plot.cLHS_result}}
#' @references *For the initial cLHS method:
#' 
#' Minasny, B. and McBratney, A.B. 2006. A conditioned Latin hypercube method
#' for sampling in the presence of ancillary information. Computers and
#' Geosciences, 32:1378-1388.
#' 
#' *For the DLHS method:
#' 
#' Minasny, B. and A. B. McBratney, A.B.. 2010. Conditioned Latin Hypercube 
#' Sampling for Calibrating Soil Sensor Data to Soil Properties. In: Proximal 
#' Soil Sensing, Progress in Soil Science, pages 111-119.
#' 
#' *For the cost-constrained implementation:
#' 
#' Roudier, P., Beaudette, D.E. and Hewitt, A.E. 2012. A conditioned Latin
#' hypercube sampling algorithm incorportaing operational constraints. In:
#' Digital Soil Assessments and Beyond. Proceedings of the 5th Golobal Workshop
#' on Digital Soil Mapping, Sydney, Australia.
#' 
#' 
#' @examples
#' 
#' df <- data.frame(
#'   a = runif(1000), 
#'   b = rnorm(1000), 
#'   c = sample(LETTERS[1:5], size = 1000, replace = TRUE)
#' )
#' 
#' # Returning the indices of the sampled points
#' res <- clhs(df, size = 50, iter = 100, progress = FALSE, simple = TRUE)
#' str(res)
#' 
#' # Returning a cLHS_result object for plotting
#' res <- clhs(df, size = 50, iter = 100, progress = FALSE, simple = FALSE)
#' str(res)
#' plot(res)
#' 
#' # Method DLHS with a linear increase of the strata weight (i.e. probability to be sampled)
#' # from 1 for the middle starta to 3 for the edge of the distribution
#' linear_increase <- 1+(2/24)*0:24
#' eta <- matrix(c(rev(linear_increase), linear_increase), ncol = 2, nrow = 50)
#' res <- clhs(df, size = 50, iter = 100, eta = eta, progress = FALSE, simple = FALSE)
#' str(res)
#' plot(res)  
#' 
#' @include clhs-data.frame.R
#' @export
clhs <- function(x, size, include,  cost,  iter,  temp,  tdecrease,  weights, eta, obj.limit, length.cycle, simple, progress, track) UseMethod("clhs")
