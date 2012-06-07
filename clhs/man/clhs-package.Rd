\name{clhs-package}
\alias{clhs-package}

\docType{package}

\title{
Conditioned Latin Hypercube Sampling
}

\description{
This package implements the conditioned Latin hypercube sampling, as published by Minasny and McBratney (2006). This method proposes to stratify sampling in presence of ancillary data.

An extension of this method, which propose to associate a cost to each individual and take it into account during the optimisation process, is also proposed (Roudier et al., 2012).
}

\author{
Pierre Roudier

Maintainer: Pierre Roudier <roudierp@landcareresearch.co.nz>
}

\references{
*For the initial cLHS method:

Minasny, B. and McBratney, A.B. 2006. A conditioned Latin hypercube method for sampling in the presence of ancillary information. Computers and Geosciences, 32:1378-1388.

*For the cost-constrained implementation:

Roudier, P., Beaudette, D.E. and Hewitt, A.E. 2012. A conditioned Latin hypercube sampling algorithm incorportaing operational constraints. In: Digital Soil Assessments and Beyond. Proceedings of the 5th Golobal Workshop on Digital Soil Mapping, Sydney, Australia.
}

\keyword{package}
\seealso{\code{sample}}

\examples{
df <- data.frame(a = runif(1000), b = rnorm(1000), c = sample(LETTERS[1:5], size = 1000, replace = TRUE))
res <- clhs(df, size = 50, iter = 2000, progress = FALSE)
str(res)
}
