#' @include clhs.R utils.R
#' @rdname clhs
#' @method clhs SpatialPointsDataFrame
#' @importFrom methods is
#' @export
#' @noRd
clhs.SpatialPointsDataFrame <- function(
  x, # data
  ...,
  use.coords = FALSE
  ){

  if (use.coords) {
    df <- data.frame(
      x@data,
      x@coords
    )
  } else {
    df <- x@data
  }
  
  spl <- clhs(x = df, ...)

  if (is(spl, "cLHS_result")) {
    spl$initial_object <- x # replacing the data.frame by the SPDF object
    spl$sampled_data <- x[spl$index_samples, ]
  }
  
  spl
}
