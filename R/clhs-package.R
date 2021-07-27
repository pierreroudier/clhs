#' Conditioned Latin Hypercube Sampling
#' 
#' This package implements the conditioned Latin hypercube sampling, as
#' published by Minasny and McBratney (2006) and the DLHS variant method 
#' (Minasny and McBratney, 2010).. This method proposes to stratify sampling in 
#' presence of ancillary data.
#' 
#' An extension of this method, which propose to associate a cost to each
#' individual and take it into account during the optimisation process, is also
#' proposed (Roudier et al., 2012).
#' 
#' @name clhs-package
#' @docType package
#' @author Pierre Roudier
#' 
#' @seealso \code{sample}
#' 
#' @references 
#' 
#' * For the initial cLHS method:
#' 
#' Minasny, B. and McBratney, A.B. 2006. A conditioned Latin hypercube method
#' for sampling in the presence of ancillary information. Computers and
#' Geosciences, 32:1378-1388.
#' 
#' *For the DLHS variant method:
#' 
#' Minasny, B. and A. B. McBratney, A.B.. 2010. Conditioned Latin Hypercube 
#' Sampling for Calibrating Soil Sensor Data to Soil Properties. In: Proximal 
#' Soil Sensing, Progress in Soil Science, pages 111-119. 
#' 
#' * For the cost-constrained implementation:
#' 
#' Roudier, P., Beaudette, D.E. and Hewitt, A.E. 2012. A conditioned Latin
#' hypercube sampling algorithm incorporating operational constraints. In:
#' Digital Soil Assessments and Beyond. Proceedings of the 5th Global Workshop
#' on Digital Soil Mapping, Sydney, Australia.
#' 
#' * For the similarity buffer prediction:
#' 
#' Brungard, C. and Johanson, J. 2015. The gate's locked! I can't get to the exact 
#' sampling spot... can I sample nearby? Pedometron, 37:8--10. 
#' 
#' @keywords sampling 
#' @import Rcpp RcppArmadillo
#' @useDynLib clhs
#' 
NULL

#' Conditioned Latin Hypercube Sampling result
#' 
#' A S3 class describing a cLHS result.
#' 
#' @name cLHS_result
#' @docType class
#' @return An object of class \code{cLHS_result} contains the following slots:
#' \item{index_samples }{a vector giving the indices of the chosen samples.}
#' \item{sampled_data }{the sampled data.frame.} \item{obj}{a vector giving the
#' evolution of the objective function throughout the Metropolis-Hastings
#' iterations.} \item{cost}{a vector giving the evolution of the cost function
#' throughout the Meropolis-Hastings iterations, if available, otherwise NULL.}
#' @author Pierre Roudier
#' @seealso \code{clhs}
NULL
