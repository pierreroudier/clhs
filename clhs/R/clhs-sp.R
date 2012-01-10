clhs.SpatialPointsDataFrame <- function(
  x, # data
  ...
  ){

  spl <- clhs(x = x@data, ...)

  if (is(spl), "cLHS_result")
    spl$initial_object <- x # replacing the data.frame by the SPDF object
  else
    res <- x[spl, ]
  res
}

setMethod("clhs", "SpatialPointsDataFrame", clhs.SpatialPointsDataFrame)
