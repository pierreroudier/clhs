clhs.Raster <- function(
  x, # data
  ...
  ){
  spdf <- as(x, "SpatialPixelsDataFrame")
  spl <- clhs(x = spdf, ...)

  spl
}

setMethod("clhs", "Raster", clhs.Raster)
