\name{clhs}
\alias{clhs}

\title{Conditioned Latin Hypercube Sampling}

\description{
Implementation of the conditioned Latin hypercube sampling, as published by Minasny and McBratney (2006). This method proposes to stratify sampling in presence of ancillary data.
}

\usage{
clhs(x, size, iter = 10000, tdecrease = 0.95, weights = list(numeric = 1, factor = 1, correlation = 1), progress = TRUE)
}

\arguments{
  \item{x}{a \code{data.frame} object.}
  \item{size}{a non-negative integer giving the number of items to choose.}
  \item{iter}{a positive number, giving the number of iterations for the Metropolis-Hastings annealing process. See `Details'.}
  \item{tdecrease}{a number betwen 0 and 1, giving the rate at which temperature decreases in the simulated annealing process. See `Details'.}
  \item{weights}{a list a length 3, giving the relative weights for continuous data, categorical data, and correlation between variables.}
  \item{progress}{TRUE or FALSE, displays a progress bar.}
}

\details{
Details about the algorithm to follow here.
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
