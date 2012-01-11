\name{clhs-package}
\alias{clhs-package}

\docType{package}

\title{
Conditioned Latin Hypercube Sampling
}

\description{
This package implements the conditioned Latin hypercube sampling, as published by Minasny and McBratney (2006). This method proposes to stratify sampling in presence of ancillary data.
}

\author{
Pierre Roudier

Maintainer: Pierre Roudier <roudierp@landcareresearch.co.nz>
}

\references{
Minasny, B. and McBratney, A.B. 2006. A conditioned Latin hypercube method for sampling in the presence of ancillary information. Computers and Geosciences, 32:1378-1388.
}

\keyword{package}
\seealso{\code{sample}}

\examples{
df <- data.frame(a = runif(1000), b = rnorm(1000), c = sample(LETTERS[1:5], size = 1000, replace = TRUE))
res <- clhs(df, size = 50, iter = 2000, progress = FALSE)
str(res)
}
