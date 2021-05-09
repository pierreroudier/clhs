#' @include clhs.R utils.R
#' @rdname clhs
#' @method clhs sf
#' @importFrom methods is 
#' @importFrom sf st_geometry_type st_coordinates st_set_geometry
#' @noRd
#' @export
clhs.sf <- function(
  x, # data
  ...,
  use.coords = FALSE
){
  
  # Check geometry type if we use coordinates
  if (use.coords) {
    
    # When coords are used only points/multipoints are supported
    if (! st_geometry_type(x, by_geometry = FALSE) %in% c("POINT", "MULTIPOINT")) {
      stop(
        "When use.coords` is set to TRUE, only POINT/MULTIPOINT geometries are supported.",
        call. = FALSE
      )
    }
    
    coords <- st_coordinates(x)
    df <- data.frame(
      st_set_geometry(x, NULL),
      coords
    )
    
  } else {
    df <- st_set_geometry(x, NULL)
  }
  
  spl <- clhs.data.frame(x = df, ...)
  
  if (is(spl, "cLHS_result")) {
    spl$initial_object <- x # replacing the data.frame by the original sf object
    spl$sampled_data <- x[spl$index_samples, ]
  }
  
  spl
}
