\name{cLHS_result}

\docType{class}

\alias{cLHS_result}

\title{Conditioned Latin Hypercube Sampling result}

\description{
A S3 class describing a cLHS result.
}

\value{
An object of class \code{cLHS_result} contains the following slots:
  \item{index_samples }{a vector giving the indices of the chosen samples.}
  \item{sampled_data }{the sampled data.frame.}
  \item{obj}{a vector giving the evolution of the objective function throughout the Meropolis-Hastings iterations.}
  \item{cost}{a vector giving the evolution of the cost function throughout the Meropolis-Hastings iterations, if available, otherwise NULL.}
}

\author{Pierre Roudier}

\seealso{\code{clhs}}
