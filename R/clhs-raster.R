#' @include clhs-data.frame.R
#' @rdname clhs
#' @importFrom raster rasterToPoints
#' @method clhs Raster
#' @export
clhs.Raster <- function(
  x, # data
  ...
  ){
  spdf <- rasterToPoints(x, spatial = TRUE)
  spl <- clhs(x = spdf, ...)
  
  dots <- list(...)
  simple <- dots$simple
  
  if (length(simple) == 0) simple <- TRUE
  
  if (!simple) {
    spl$initial_object <- x
  }
  
  spl
}
