\name{plot.cLHS_result}
\alias{plot.cLHS_result}
\title{Plot cLHS results}
\usage{\method{plot}{cLHS_result}(x, modes = "obj", ...)
}

\arguments{
\item{x}{Object of class \dQuote{cLHS_result}.}
\item{modes}{A character vector describing the plot to produce (see Details)}
\item{...}{Other ggplot2 plotting parameters.}
}

\description{Produces a plot illustrating the result of a cLHS sampling procedure. }

\value{A ggplot object}

\details{
The subplots to be included in the final illustration are controlled by the \code{mdoe} option:
- \code{"obj"} adds the evolution of the objective function over the iterations
- \code{"hist"} adds the comparison of the distributions of each variables in both the original object and the sampled result.
}

\seealso{\code{\link{clhs}}}

\author{Pierre Roudier}

\examples{
df <- data.frame(a = runif(1000), b = rnorm(1000), c = sample(LETTERS[1:5], size = 1000, replace = TRUE))
res <- clhs(df, size = 50, iter = 2000, progress = FALSE, simple = FALSE)
# You can plot only the objective function
plot(res, mode = "obj")
# Or you can compare the distribution in the original object 
# and in the sampled result
plot(res, mode = c("obj", "hist"))
}
