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
  \item{obj_function }{a vector giving the evolution of the objective function throughout the Meropolis-Hastings iterations.}
}

\author{Pierre Roudier}

\seealso{\code{clhs}}
