clhs.Raster <- function(
  x, # data
  ...
  ){
  spdf <- rasterToPoints(x, spatial = TRUE)
  spl <- clhs(x = spdf, ...)
  
  dots <- list(...)
  simple <- dots$simple
  
  if (!simple) {
    spl$initial_object <- x
  }
  
  spl
}

setMethod("clhs", "Raster", clhs.Raster)
