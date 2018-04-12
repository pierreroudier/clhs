## TODO: convert to gower package: https://cran.r-project.org/web/packages/gower/

#' @name Gower similarity
#' @description Calculates Gower's similarity index for every pixel within an given radius buffer of each sampling point
#' @param covs raster stack of environmental covariates
#' @param pts sampling points, object of class SpatialPointsDataframe
#' @param buffer Radius of the disk around each point that similarity will be calculated
#' @param fac numeric, can be > 1, (e.g., fac = c(2,3)). Raster layer(s) which are categorical variables. Set to NA if no factor is present
#' @param ... passed to plyr::alply
#' @return a RasterStack
#' @author Colby Brungard (cbrung@nmsu.edu)
#' 
#' @importFrom stats complete.cases
#' @importFrom raster extract cellFromXY raster stack nlayers
#' @importFrom cluster daisy
#' @importFrom plyr alply
#' 
#' @export
#' 
#' @example 
#' library(raster)
#' library(sp)
#' 
#' data(meuse.grid)
#' coordinates(meuse.grid) = ~x+y
#' proj4string(meuse.grid) <- CRS("+init=epsg:28992")
#' gridded(meuse.grid) = TRUE
#' ms <- stack(meuse.grid)
#' set.seed(1)
#' pts <- clhs(ms, size = 3, iter = 100, progress = FALSE, simple = TRUE)
#' gw <- similarity_buffer(ms, pts, buffer = 500)
#' plot(gw)
#' 
similarity_buffer <- function(covs, pts, buffer, fac = NA, ...) {

  # Iterate over every point
  res_l <- plyr::alply(pts, 1, function(coords) {
    
    # 2. Extract all cells within x m of the sampling points. 
    buff_data <- extract(
      x = covs, 
      y = coords, 
      buffer = buffer, 
      cellnumbers = TRUE, 
      method = 'simple', 
      df = TRUE
    )
    
    # 3. Apply gowers similarity index to each element of list of extracted raster values
    
    # Get the cell numbers from each sample point to identity the right column in the similarity matrix. 
    cellnum <- cellFromXY(covs, coords)
    
    # 3.b Calculate Gower's similarity index around each point. 
    #   I used Gower's because it can handle categorical covariates, 
    #   but there could be other options. 
    
    # Only retain cases without NA values
    buff_data <- data.frame(buff_data[complete.cases(buff_data), ])
    
    # If there are some factor data
    if (!any(is.na(fac))) {
      buff_data[, fac + 1] <- lapply(buff_data[fac + 1], factor)
    }

    # Calculate gowers similarity index 
    gower_dissim <- daisy(x = buff_data[, names(covs)], metric = 'gower', warnBin = FALSE)
    # turn dissimilarity object to matrix
    gower_dissim <- cbind(buff_data$cell, as.matrix(gower_dissim)) 
    
    # Select the row of similarity indices with cell number equal to the cell number of the 
    # sample point and convert dissimilarity to similarity by subtracting from 1.  
    gower_sim <- 1 - gower_dissim[gower_dissim[, 1] == cellnum, ] 
    
    # Combine the cellnumbers of the raster to the similarity index. 
    res_df <- data.frame(cellnum = buff_data$cell, similarity = gower_sim[-1])
    
    # Define the raster layer storing the results
    res_r <- raster(covs, layer = 0) 
    
    # Index the original raster by the cell numbers and replace the NA values with the similarity values. 
    # This results in a raster with similarity values in the buffers around each point and NA everywhere else. 
    res_r[res_df$cellnum] <- res_df$similarity
    
    res_r
  }, ...)
  
  res_s <- stack(res_l)
  names(res_s) <- paste0('similarity_point_', 1:nlayers(res_s))
  
  res_s
} 