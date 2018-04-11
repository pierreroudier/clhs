# This code written by Colby Brungard to produce Gowers similarity indecies around sampling points. 
# 11/13/2017

# Colby Brungard, PhD
# Assistant Professor of Pedology | Plant and Environmental Sciences Dept.  
# New Mexico State University | Las Cruces, NM 88003
# +1-575-646-1907
# cbrung@nmsu.edu

## TODO: convert to gower package: https://cran.r-project.org/web/packages/gower/

## TODO: ouput vector of file names for later use


# Function to calculate Gower's similarity index for every pixel within an X m radius buffer of each sampling point
# 

#' @name Gower similarity
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
  # rstack: raster stack of environmental covariates
  # samps: SpatialPointsDataframe of your sampling points. 
  # This MUST have one column named 'ID' of sample point names to name the output rasters
  # Buffer: numeric. Buffer distance around each point that similarity will be calculated. 
  # If buffer = 250 the final buffer will be 500 m in diameter.
  # fac: numeric, can be > 1, (e.g., fac = c(2,3)). Raster layer(s) which are categorical variables. Set to NA if no factor is present. 
  
  # Iterate over every point. This keeps memory usage small
  f.list <- list() #vector(mode = 'character', length = nrow(samps))
  
  for(i in 1:nrow(samps)){
    
    #2. Extract all cells within x m of the sampling points. 
    buff <- extract(x = rstack, y = samps[i,], buffer = buffer, cellnumbers = TRUE, method = 'simple')
    
    # 3.Apply gowers similarity index to each element of list of extracted raster values.
    # Get the cell numbers from each sample point to identity the right column in the similarity matrix. 
    cellNum <- extract(x = rstack, y = samps[i,], method = 'simple', cellnumbers = TRUE, small = TRUE)
    cellnum <- cellNum[,1] #just need the cell numbers. 
    
    # 3.b Function to calculate Gower's similarity index around each point. I used Gower's because it can handle categorical covariates, but there are other options. 
    
    # First convert each list element in buff to a data frame, then remove NA values. If values are on edges this makes the elements of similarL not the same so they can't be rbind-ed together, third Force any factors to be a factor (+1 because the cell numbers are added as the first column)
    tL <- data.frame(buff) 
    tL <- tL[complete.cases(tL),]
    
    # Warnings are suppressed because if fac is > 1 the is.na function warns that only the first value will be used.
    # This does not affect the output whatsoever.   
    # All other elements of the result list are as expected.
    suppressWarnings(
      if (is.na(fac)) {
        tL <- tL
      } else {
        tL[,fac+1] <- lapply(tL[fac + 1], factor)
      })
    
    # Calculate gowers similarity index and make dissimilarity object a matrix
    GL <- daisy(x = tL[, -1], metric = 'gower')
    
    gmatL <- cbind(tL$cell, as.matrix(GL)) 
    
    # Select the row of similarity indices with cell number equal to the cell number of the sample point and convert dissimilarity to similarity by subtracting from 1.  
    finalgL <- 1 - gmatL[gmatL[, 1] == cellnum,] 
    
    # Combine the cellnumbers of the raster to the similarity index. 
    similarL <- list()
    similarL[[1]] <- cbind(tL$cell, finalgL[-1]) 
    
    # Note this may return several warning messages from the daisy package. See https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/daisy.html for more information on these warnings.
    
    # Convert list of cell numbers and similarity values to a dataframe
    s.df <- do.call(rbind.data.frame, similarL) 
    names(s.df) <- c('CellNum', 'Similar') #Give better names
    
    # Define a raster to hold the similarity values by selecting the first layer of the raster stack, and make all cells NA
    r <- raster(rstack, layer = 0) 
    
    # Index the original raster by the cell numbers and replace the NA values with the similarity values. 
    # This results in a raster with similarity values in the buffers around each point and NA everywhere else. 
    r[s.df$CellNum] <- s.df$Similar
    names(r) <- paste0('SimilarityIndex_', i)
    
    f.list[[i]] <- r
  }
  
  stack(f.list)
} 