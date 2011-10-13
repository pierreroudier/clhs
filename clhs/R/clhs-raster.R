clhs.Raster <- function(
  x, # data
  ...
  ){
  spdf <- as(x, "SpatialPixelsDataFrame")
  spl <- clhs(x = spdf, ...)

  if(ncol(spl) > 1)
    res <- stack(spl)
  else
    res <- raster(spl)

  res
}

setMethod("clhs", "Raster", clhs.Raster)
