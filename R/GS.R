## TODO: convert to gower package: https://cran.r-project.org/web/packages/gower/

# 
# 

#' @name Gower similarity
#' @description Calculates Gower's similarity index for every pixel within an given radius buffer of each sampling point
#' @param rstack raster stack of environmental covariates
#' @param samps sampling points, object of class SpatialPointsDataframe
#' @param buffer Radius of the disk around each point that similarity will be calculated
#' @param fac numeric, can be > 1, (e.g., fac = c(2,3)). Raster layer(s) which are categorical variables. Set to NA if no factor is present
#' @return a RasterStack
#' @author Colby Brungard (cbrung@nmsu.edu)
#' @example 
#' library(raster)
#' library(sp)
#' 
#' slogo <- stack(system.file("external/rlogo.grd", package="raster")) 
#' coords <- data.frame(x = sample(1:101, size = 25), y = sample(1:77, size = 25))
#' spdf <- SpatialPointsDataFrame(coords, data = data.frame(ID = 1:25))
#' 
#' gw <- GS(rstack = slogo, samps = spdf, buffer = 25, fac = NA)
#' plot(gw)
#' 
GS <- function(rstack, samps, buffer, fac = NA) {

  # Iterate over every point
  # res <- plyr::alply(coordinates(samps), function(coords) {
  #   
  #   # Extract all cells within x m of the sampling points. 
  #   buff <- extract(x = rstack, y = coords, buffer = buffer, cellnumbers = TRUE, method = 'simple', df  =TRUE)
  #   
  #   # Apply Gower's similarity to each element of the list of extracted raster values.
  #   # Get the cell numbers from each sample point to identity the right column in the similarity matrix. 
  #   cellNum <- extract(x = rstack, y = coords, method = 'simple', cellnumbers = TRUE, small = TRUE)
  #   cellnum <- cellNum[, 1] #just need the cell numbers. 
  #   
  # })
  
  # Iterate over every point. This keeps memory usage small
  f.list <- list()
  
  for(i in 1:nrow(samps)){
    
    # 2. Extract all cells within x m of the sampling points. 
    buff_data <- extract(x = rstack, y = samps[i, ], buffer = buffer, cellnumbers = TRUE, method = 'simple', df = TRUE)
    
    # 3. Apply gowers similarity index to each element of list of extracted raster values
    
    # Get the cell numbers from each sample point to identity the right column in the similarity matrix. 
    cellnum <- cellFromXY(rstack, samps[i,])
    
    # 3.b Calculate Gower's similarity index around each point. 
    #   I used Gower's because it can handle categorical covariates, 
    #   but there could be other options. 
    
    # Only retain cases without NA values
    buff_data <- data.frame(buff_data[complete.cases(buff_data), ])
    
    # Warnings are suppressed because if fac is > 1 the is.na function 
    # warns that only the first value will be used.
    # This does not affect the output whatsoever.   
    # All other elements of the result list are as expected.
    suppressWarnings(
      if (!is.na(fac)) {
        buff_data[, fac + 1] <- lapply(buff_data[fac + 1], factor)
      }
    )
    
    # Calculate gowers similarity index and make dissimilarity object a matrix
    gower_dissim <- daisy(x = buff_data[, names(rstack)], metric = 'gower')
    gower_dissim <- cbind(buff_data$cell, as.matrix(gower_dissim)) 
    
    # Select the row of similarity indices with cell number equal to the cell number of the 
    # sample point and convert dissimilarity to similarity by subtracting from 1.  
    gower_sim <- 1 - gower_dissim[gower_dissim[, 1] == cellnum, ] 
    
    # Combine the cellnumbers of the raster to the similarity index. 
    s.df <- data.frame(cellnum = buff_data$cell, similarity = gower_sim[-1])
    
    # Define the raster layer storing the results
    r <- raster(rstack, layer = 0) 
    
    # Index the original raster by the cell numbers and replace the NA values with the similarity values. 
    # This results in a raster with similarity values in the buffers around each point and NA everywhere else. 
    r[s.df$cellnum] <- s.df$similarity
    
    names(r) <- paste0('SimilarityIndex_', i)
    
    f.list[[i]] <- r
  }
  
  stack(f.list)
} 