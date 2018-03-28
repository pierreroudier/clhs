\name{GS}
\alias{GS}

\title{GS}
\description{GS}
\usage{
GS(rstack, samps, buffer, fac)
}

\arguments{
  \item{rstack}{A raster stack of environmental covariates, projected CRS.}
  \item{samps}{SpatialPointsDataframe of your sampling points, projected CRS, MUST have one column named 'ID' of sample point names to name the output rasters.}
  \item{buffer}{Buffer distance around each point that similarity will be calculated, units of CRS.}
  \item{fac}{numeric, can be > 1, (e.g., fac = c(2,3)). Raster layer(s) which are categorical variables. Set to NA if no factor is present.}
}

\details{The size of the buffer will have an effect on the calculated similarity values. This is because Gower's similariy index standardizes each variable (all similarity algorithms do*) by dividing by the range of the values. Selecting a larger buffer size will result in a larger range of values thus resulting in different similarity values. A larger buffer will result in higher similarity values at any given pixel. In my experience this becomes less significant when a buffer is > 150 m... One solution to this would be to again standardize each similarity raster to the range [0,1] using the range of the largest buffer, but I don't think that it is that important because the pattern of which cells have the highest values remains the same, while only the actual similarity values change. It is probably best to put some thought into what buffer size is most appropriate.}

\value{...}

\references{...}

\author{Colby Brungard}

\note{If a raster is poorly symbolized when viewed in a GIS (e.g., the lower range is listed as nan), manually rescale the symbology from 0.5 to 1. The algorithm worked, it is just the GIS that is having problems. }


\examples{
\dontrun{
library(raster)
library(rgdal)
library(cluster)

# Define input variables
rs <- stack("./rasterstack.tif")
s <- shapefile("./Samples.shp")

# Run the function. This takes about 2 seconds per sample on a 2.6 GHz CPU
GS(rstack = rs, samps = s, buffer = 250, fac = c(2,3))

}
}

