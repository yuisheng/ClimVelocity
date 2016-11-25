library(rgdal)
library(raster)
library(sp)
library(proj4)
library(rasterVis)
library(doParallel)
library(foreach)

setwd("D:/Eawag/course - Spatial Modelling/ClimateVelocity/")

# ch <- read.csv('inputs/Swiss_border_coordinates.dat', sep='\t', header=TRUE, stringsAsFactors = FALSE)
# ch <- ch[, c("POINT_X", "POINT_Y")]
# Current climate must be divided by 100 to get degrees Celsius
cc <- stack("inputs/current_climate/currentclimate_1971-2000.grd")
cc <- crop(cc, extent(c(5e5, 1e6, 0, 5e5)))

fc1 <- stack("inputs/future_climate/rcp45/CLMcom_CCLM4-8-17_MOHC_HadGEM2-ES_ar5_wc_rcp45.grd")
# fc2 <- stack("inputs/future_climate/rcp45/DMI_HIRHAM5_ICHEC_EC-EARTH_ar5_wc_rcp45.grd")
# fc3 <- stack("inputs/future_climate/rcp45/KNMI_RACMO22E_ICHEC_EC-EARTH_ar5_wc_rcp45.grd")

fc1 <- crop(fc1, extent(c(5e5, 1e6, 0, 5e5)))
# fc2 <- crop(fc2, extent(c(5e5, 1e6, 0, 5e5)))
# fc3 <- crop(fc3, extent(c(5e5, 1e6, 0, 5e5)))

### FUNCTION: ToCell
# Purpose: ToCell() takes one present and one future climate raster, one climate variable of interest, and a desired climate change "window" or "accuracy" to compute the future climate 
# cell with the nearest value to every present climate cell. 

# Implementation: decomposes rasters to matrices and matrices to vectors, 
#   uses indices, and loops through present climate vectors in parallel 

# Arguments:
# - pRaster: cropped present climate raster
# - fRaster: cropped future climate raster
# - climVar: string of climate variable: tas or prcp
# - accuracy: future cells with value within present cell value +/- accuracy (integer)
# - cores: number of physical CPU cores to be used in computation (integer)

# Outputs:
# Outputs list of two matrix objects of equal dimensions:
# - "fromcell" contains x/y/climVar of current climate cells (vector start point)
# - "tocell" contains x/y/climVar of future cells corresponding to current climate cells (vector endpoint)
ToCell <- function(pRaster, fRaster, climVar, accuracy, cores) {
  # Convert input rasters to matrices
  p <- rasterToPoints(pRaster)
  f <- rasterToPoints(fRaster)
  
  # Correct temperature and precipitation values
  p[, "tas"] <- round((p[, "tas"]/100), digits = 1)
  f[, "prcp"] <- round((f[, "prcp"]*10), digits = 1)
  
  # Decompose input matrices into vectors
  # Target climate vectors
  upper <- as.vector(p[, climVar] + accuracy)
  lower <- as.vector(p[, climVar]) - accuracy
  # Future climate vectors
  fx <- as.vector(f[, "x"])
  fy <- as.vector(f[, "y"])
  fvar <- as.vector(f[, climVar])
  
  # Prepare parallelized loop
  # Create cluster of worker processes
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  # Determine number of iterations
  l <- nrow(p)
  
  # Distribute loop to workers, rbind results
  out <- foreach(i = 1:l, .combine=rbind) %dopar% {
    # get the median future climate value within temperature limits
    target <- median(fvar[fvar > lower[i] & fvar < upper[i]])
    
    if (is.na(target) == TRUE){
      to.cell <- c(NA, NA, NA)
    }
    else {
      # within the future temperature, locate the cell closest to the target
      index <- which.min(abs(fvar - target))
      rm(target)
      
      # output the future cell information (x, y, climVar)
      to.cell <- c(fx[index], fy[index], fvar[index])
    }
  }
  stopCluster(cl) # Stop the cluster
  
  rownames(out) <- 1:nrow(out)
  colnames(out) <- c("x", "y", climVar)
  
  # Output list of matrices: fromcell for current climate, tocell for climate destination  
  result <- list("fromcell" = p[, c("x", "y", climVar)], "tocell" = out)
  return(result)
}

x <- system.time(ToCell(cc, fc1, "tas", 1, 4))

colnames(x$fromcell) <- c("x1", "y1", "tas1")
colnames(x$tocell) <- c("x2", "y2", "tas2")
sp <- cbind(x$fromcell, x$tocell)


f <- rasterToPoints(fc1)
# Correct temperature and precipitation values
f[, "prcp"] <- f[, "prcp"]*10


# # Plot first 10 cells using XY limits
# plot(x[1,"x1"], x[1,"x1"], xlim=c(min(x[,"x1"]), max(x[,"x1"])), ylim = c(min(x[,"y1"]), max(x[,"y1"])))
# 
# # Take the first ten cells and plot their vectors
# i <- x[1:10,]
# arrows(i[,1], i[,2], i[,4], i[,5])

















