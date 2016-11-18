library(rgdal)
library(raster)
library(sp)
library(proj4)
library(doParallel)
library(foreach)

setwd("D:/Eawag/course - Spatial Modelling/ClimateVelocity/ClimVelocity")

# Current climate must be divided by 100 to get degrees Celsius
cc <- stack("current_climate/currentclimate_1971-2000.grd")

fc1 <- stack("future_climate/rcp45/CLMcom_CCLM4-8-17_MOHC_HadGEM2-ES_ar5_wc_rcp45.grd")
fc2 <- stack("future_climate/rcp45/DMI_HIRHAM5_ICHEC_EC-EARTH_ar5_wc_rcp45.grd")
fc3 <- stack("future_climate/rcp45/KNMI_RACMO22E_ICHEC_EC-EARTH_ar5_wc_rcp45.grd")


# present <- rasterToPoints(cc)
# present[, "tas"] <- present[, "tas"]/100 # Get degrees Celsius
# f <- rasterToPoints(fc1)


# For every current climate cell, get future climate cells within climate window
# Apply() approach:
# upper/lower temp values calculated within each iteration
# all objects are deleted ASAP
#
# p <- as.vector(present[,3])
# 
# sapply(p, function(x){
#   # define temperature limits
#   upper <- x + accuracy
#   lower <- x - accuracy
#   
#   # get median future climate that falls within temperature limits of cell
#   target <- median(f[ f[ ,3] < upper & f[ ,3] > lower, 1:3])
#   
#   rm(upper, lower)
#   # within the future temperature, get the cell closest to the target
#   to.cell <- f[which.min(abs(f[, 3] - target)), 1:3]
#   
#   rm(target)
#   return(to.cell)
# })

# ### LOOP approach ###:
# Upper/lower temp values calculated beforehand, added to present climate matrix
# Column names never used; index locations used for performance
# Length of vector calculated only once
# 
# calculate upper and lower temperature per cell
present <- rasterToPoints(cc)
present[, "tas"] <- present[, "tas"]/100 # Get degrees Celsius
f <- rasterToPoints(fc1)

accuracy <- 1
px <- as.vector(present[, 1])
py <- as.vector(present[, 2])
ptemp <- as.vector(present[, 3])
pprecip <- as.vector(present[, 4])
upper <- ptemp + accuracy
lower <- ptemp - accuracy
# plist <- list(px, py, ptemp, pprecip, upper, lower) 
# plist <- lapply(seq_len(ncol(present)), function(i) present[,i])

fx <- as.vector(f[, 1])
fy <- as.vector(f[, 2])
ftemp <- as.vector(f[, 3])
fprecip <- as.vector(f[, 4])

# Pre-allocate matrix containing results
# d <- matrix(nrow=nrow(present), ncol=3) # preallocate matrix

l <- length(ptemp) # 2,695,527 elements
for (i in 1:l){
  # get future climate that falls within temperature limits of cell
  cat(i, "\n")
  target <- median(ftemp[ftemp > lower[i] & ftemp < upper[i]])
  if (is.na(target) == TRUE){
    to.cell <- c(NA, NA, NA)
  }
  
  else {
    # within the future temperature, locate the cell closest to the target
    system.time(index <- which.min(abs(ftemp - target)))
    
    rm(target)
    # Retrieve the information for the located cell
    
    to.cell <- c(fx[index], fy[index], ftemp[index])
  }
}

### Parallelized loop ###
cl <- makeCluster(4)
registerDoParallel(cl)

l <- length(ptemp) # 2,695,527 elements
system.time(output <- foreach(i = 1:l, .combine=rbind) %dopar% {
  # get the median future climate value within temperature limits
  target <- median(ftemp[ftemp > lower[i] & ftemp < upper[i]])
  if (is.na(target) == TRUE){
    to.cell <- c(NA, NA, NA)
  }
  else {
    # within the future temperature, locate the cell closest to the target
    index <- which.min(abs(ftemp - target))
    rm(target)
    
    # output the future cell information
    to.cell <- c(fx[index], fy[index], ftemp[index])
  }
})
rownames(output) <- 1:nrow(output)
stopCluster(cl)

### Compiled function ###
library(compiler)

CV <- function(ptemp, ftemp, lower, upper){
  l <- length(ptemp) # 2,695,527 elements
  for (i in 1:l){
    # get future climate that falls within temperature limits of cell
    target <- median(ftemp[ftemp > lower[i] & ftemp < upper[i]])
    if (is.na(target) == TRUE){
      to.cell <- c(NA, NA, NA)
    }
    else {
      # within the future temperature, locate the cell closest to the target
      index <- which.min(abs(ftemp - target))
      rm(target)
      # Retrieve the information for the located cell
      cat(i, "\n")
      to.cell <- c(fx[index], fy[index], ftemp[index])
    }
  }
}

CVc <- cmpfun(CV)



### SCRATCH ###
# For cell x in rasterA
# define accuracy window
# get distances of cells in rasterB within accuracy window of cell x
# take the mean / median distance to choose a cell in rasterB
# calculate the vector	

# <21781> +title=CH1903 / LV03 +proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs

# # WGS 1984
# crs.new <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# cc.new <- projectRaster(cc, crs = crs.new)
# # Great sphere distance between lat/long points with cosine
# distm(c(6.376007, 52.33323), c(6.388707, 52.33323), fun=distCosine)
# 
# # Ellipsoidal distance along lat/long points with Bessel 1841
# distGeo(c(6.376007, 52.33323), c(6.388707, 52.33323), a = 6377397.155, f = 1/299.1528434)
# 
# test <- spTransform(window[,1:2], crs.new)


