library(rgdal)
library(raster)
library(sp)
library(proj4)
library(rasterVis)
library(doParallel)
library(foreach)
library(dplyr)
library(SDMTools) 
setwd("D:/Eawag/course - Spatial Modelling/ClimateVelocity/")

# Current climate must be divided by 100 to get degrees Celsius
cc <- stack("inputs/current_climate/currentclimate_1971-2000.grd")

# Scenario 1: RCP45
rcp45.fc1 <- stack("inputs/future_climate/rcp45/CLMcom_CCLM4-8-17_MOHC_HadGEM2-ES_ar5_wc_rcp45.grd")
rcp45.fc2 <- stack("inputs/future_climate/rcp45/DMI_HIRHAM5_ICHEC_EC-EARTH_ar5_wc_rcp45.grd")
rcp45.fc3 <- stack("inputs/future_climate/rcp45/KNMI_RACMO22E_ICHEC_EC-EARTH_ar5_wc_rcp45.grd")

# Scenario 2: RCP85
rcp85.fc1 <- stack("inputs/future_climate/rcp85/CLMcom_CCLM4-8-17_MOHC_HadGEM2-ES_ar5_wc_rcp85.grd")
rcp85.fc2 <- stack("inputs/future_climate/rcp85/DMI_HIRHAM5_ICHEC_EC-EARTH_ar5_wc_rcp85.grd")
rcp85.fc3 <- stack("inputs/future_climate/rcp85/KNMI_RACMO22E_ICHEC_EC-EARTH_ar5_wc_rcp85.grd")

### FUNCTION ###
# Climate Velocity (CV)
# Arguments:
# - pRaster: present climate raster
# - fRaster: future climate raster
# - climVar: string of current climate variable: "tas" / "prcp"
# - threshold: future cells with value within present cell value +/- threshold (integer)
# - cores: number of physical CPU cores to be used in computation (integer)

CV <- function(pRaster, fRaster, cVar, cores, threshold) {
  present <- rasterToPoints(pRaster)
  present[,"tas"] <- present[,"tas"]/100
  present <- as.data.frame(present)
  
  future <- rasterToPoints(fRaster)
  future[,"prcp"] <- future[,"prcp"]*10
  future <- as.data.frame(future)
  colnames(future) <- c("x", "y", "ftas", "fprcp")
  
  # Set future climate variable name dependent on cVar
  if (cVar == "tas"){
    fcVar <- "ftas"
  }
  if(cVar == "prcp"){
    fcVar <- "fprcp"
  }
  
  ### SPATIAL GRADIENT ###
  t <- threshold          # plus/minus threshold to define climate match
  t <- 1/(t*2)            # inverse for rounding, double for plus/minus
  
  x <- present$x                    # vector of grid cell x coordinates
  y <- present$y                    # vector of grid cell y coordinates
  p <- round(present[, cVar]*t)/t   # vector of rounded present climate values 
  f <- round(future[, fcVar]*t)/t   # vector of rounded future climate values
  
  u <- unique(p)[order(unique(p))]            # list of unique climate values in p
  m <- sapply(u, function(i){c(which(i==f))}) # list of climate matches for unique values
  
  ### PARALLELIZED LOOP ###
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  system.time(d <- foreach(i = 1:length(p), .combine='c') %dopar% {
    mi   <- m[[which(u==p[i])]]          # recalls list of climate matches for p[i]
    sqrt(min((x[i]-x[mi])^2 + (y[i]-y[mi])^2))    # distance to closest match
  })
  
  stopCluster(cl) # Stop the cluster

  # Keep original results and replace Inf with 10,000m
  q <- d
  q[q==Inf] <- 10000
  q <- q/1000/100 # Convert from m/100yr to km/yr
  
  df <- data.frame(x, y, round(present[,cVar]), q); colnames(df) <- c("x", "y", cVar, "dist")
  df[df$dist == 0, ] <- 0.01
  df <- left_join(df, future, by = c("x", "y")) # Join future with extra cell to present data to obtain future climvar
  df <- na.omit(df)
  
  ### TEMPORAL GRADIENT ###
  # future - current tas
  # divided by 100
  df$tempg <- (df[, fcVar] - df$tas)/100 # Divide by years 2100-2000 
  
  ### SPATIAL GRADIENT ###
  # Get the spatial gradient (degrees C/km)
  df$spatg <- ((df[, fcVar] - df$tas)*100)/df$dist
  
  # Get the temporal gradient (degrees C/yr)
  df$vector <- df$tempg/df$spatg
   
  ### PLOT VECTOR FIELD MAP ###
  proj <- CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000
              +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")
  vdf <- select(df, x, y, vector)
  r <- rasterFromXYZ(vdf, crs=proj)
  vectorplot(r, par.settings=RdBuTheme())
  
  output <- list(d, df)
  return(output)
}

# Temperature
# RCP45
t451 <- CV(cc, rcp45.fc1, "tas", 4, 0.25)
t452 <- CV(cc, rcp45.fc2, "tas", 4, 0.25)
t453 <- CV(cc, rcp45.fc3, "tas", 4, 0.25)
# RCP85
t851 <- CV(cc, rcp85.fc1, "tas", 4, 0.25)
t852 <- CV(cc, rcp85.fc2, "tas", 4, 0.25)
t853 <- CV(cc, rcp85.fc3, "tas", 4, 0.25)

# Precipitation
# RCP45
p451 <- CV(cc, rcp45.fc1, "prcp", 4, 5)
p452 <- CV(cc, rcp45.fc2, "prcp", 4, 5)
p453 <- CV(cc, rcp45.fc3, "prcp", 4, 5)
# RCP85
p851 <- CV(cc, rcp85.fc1, "prcp", 4, 5)
p852 <- CV(cc, rcp85.fc2, "prcp", 4, 5)
p853 <- CV(cc, rcp85.fc3, "prcp", 4, 5)




