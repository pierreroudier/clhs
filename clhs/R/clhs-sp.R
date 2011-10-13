clhs.SpatialPointsDataFrame <- function(
  x, # data
  ...
  ){

  spl <- clhs(x = x@data, ...)

  x[spl$index_samples, ]

}

setMethod("clhs", "SpatialPointsDataFrame", clhs.SpatialPointsDataFrame)
