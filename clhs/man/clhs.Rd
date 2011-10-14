\name{clhs}

\docType{methods}

\alias{clhs}
\alias{clhs.data.frame}
\alias{clhs.SpatialPointsDataFrame}
\alias{clhs.Raster}
\alias{clhs,data.frame-method}
\alias{clhs,Raster-method}
\alias{clhs,SpatialPointsDataFrame-method}

\title{Conditioned Latin Hypercube Sampling}

\description{
Implementation of the conditioned Latin hypercube sampling, as published by Minasny and McBratney (2006). This method proposes to stratify sampling in presence of ancillary data.
}

\usage{
clhs(x, ...)
}

\arguments{
  \item{x}{A \code{data.frame} object.}
  \item{...}{Additional arguments, see under Details.}
}

\details{
Below are the additional arguments that can be used:

- \code{size}: A non-negative integer giving the number of items to choose.

- \code{iter}: A positive number, giving the number of iterations for the Metropolis-Hastings annealing process.

- \code{tdecrease}: A number betwen 0 and 1, giving the rate at which temperature decreases in the simulated annealing process.

- \code{weights}: A list a length 3, giving the relative weights for continuous data, categorical data, and correlation between variables.

- \code{progress}: TRUE or FALSE, displays a progress bar.

}

\value{
An object of class \code{cLHS_result}, with the following elements:
\item{index_samples }{a vector giving the indices of the chosen samples.}
\item{sampled_data }{the sampled data.frame.}
\item{obj_function }{a vector giving the evolution of the objective function throughout the Meropolis-Hastings iterations.}
}

\references{
Minasny, B. and McBratney, A.B. 2006. A conditioned Latin hypercube method for sampling in the presence of ancillary information. Computers and Geosciences, 32:1378-1388.
}

\author{Pierre Roudier}

\note{}

\seealso{}

\examples{
df <- data.frame(a = runif(1000), b = rnorm(1000), c = sample(LETTERS[1:5], size = 1000, replace = TRUE))
res <- clhs(df, size = 50, iter = 2000, progress = FALSE)
str(res)
}
