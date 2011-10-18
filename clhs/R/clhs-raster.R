clhs.Raster <- function(
  x, # data
  ...
  ){
  spdf <- rasterToPoints(x, spatial = TRUE)
  spl <- clhs(x = spdf, ...)

  spl

}

setMethod("clhs", "Raster", clhs.Raster)
