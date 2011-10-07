\name{clhs-package}
\alias{clhs-package}
\alias{clhs}

\docType{package}

\title{
Conditioned Latin Hypercube Sampling
}

\description{
This package implements the conditioned Latin hypercube sampling, as published by Minasny and McBratney (2006). This method proposes to stratify sampling in presence of ancillary data.
}

\details{
\tabular{ll}{
Package: \tab clhs\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-10-07\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}

~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}

\author{
Pierre Roudier

Maintainer: Pierre Roudier <roudierp@landcareresearch.co.nz>
}

\references{
Minasny, B. and McBratney, A.B. 2006. A conditioned Latin hypercube method for sampling in the presence of ancillary information. Computers and Geosciences, 32:1378-1388.
}

~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
data(mtcars)
res <- clhs(mtcars, size = 5, iter = 100)
plot(res)
str(res)
}
