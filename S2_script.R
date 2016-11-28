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
# cc <- crop(cc, extent(c(5e5, 1e6, 0, 5e5)))

fc1 <- stack("inputs/future_climate/rcp45/CLMcom_CCLM4-8-17_MOHC_HadGEM2-ES_ar5_wc_rcp45.grd")
# fc2 <- stack("inputs/future_climate/rcp45/DMI_HIRHAM5_ICHEC_EC-EARTH_ar5_wc_rcp45.grd")
# fc3 <- stack("inputs/future_climate/rcp45/KNMI_RACMO22E_ICHEC_EC-EARTH_ar5_wc_rcp45.grd")

# fc1 <- crop(fc1, extent(c(5e5, 1e6, 0, 5e5)))

# present <- raster::mask(cc, fc1)
present <- rasterToPoints(cc)
present[,"tas"] <- present[,"tas"]/100
present <- as.data.frame(present)

future <- rasterToPoints(fc1)
future[,"prcp"] <- future[,"prcp"]*10
future <- as.data.frame(future)
colnames(future) <- c("x", "y", "ftas", "fprcp")

### SPATIAL GRADIENT ###
t <- 0.25               # plus/minus threshold to define climate match
t <- 1/(t*2)            # inverse for rounding, double for plus/minus

x <- present$x                    # vector of grid cell x coordinates
y <- present$y                    # vector of grid cell y coordinates
p <- round(present$tas*t)/t     # vector of rounded present climate values 
f <- round(future$ftas*t)/t     # vector of rounded future climate values

u     <- unique(p)[order(unique(p))]    # list of unique climate values in p
m     <- sapply(u, function(i){c(which(i==f))}) # list of climate matches for unique values

### PARALLELIZED LOOP ###
cl <- makeCluster(4)
registerDoParallel(cl)
system.time(d <- foreach(i = 1:length(p), .combine='c') %dopar% {
  mi   <- m[[which(u==p[i])]]          # recalls list of climate matches for p[i]
  sqrt(min((x[i]-x[mi])^2 + (y[i]-y[mi])^2))    # distance to closest match
})
stopCluster(cl) # Stop the cluster

# ### S2 LOOP ###
# for(i in 1:length(p)){                  # loop for all grid cells of p
#   mi   <- m[[which(u==p[i])]]          # recalls list of climate matches for p[i]
#   d[i] <- sqrt(min((x[i]-x[mi])^2 + (y[i]-y[mi])^2))    # distance to closest match
# }


# Keep original results and replace Inf with 10,000m
q <- d
q[q==Inf] <- 10000
q <- q/1000/100 # Convert to meters to km/year

df <- data.frame(x, y, round(present$tas), q); colnames(df) <- c("x", "y", "tas", "dist")
df[df$dist == 0, ] <- 0.01
df <- left_join(df, future, by = c("x", "y")) # Join future with extra cell to present data to obtain future climvar
df <- na.omit(df)

### TEMPORAL GRADIENT ###
# future - current tas
# divided by 100
df$tempg <- (df$ftas - df$tas)/100 # Divide by years 2100-2000 

### SPATIAL GRADIENT ###
# Get the spatial gradient (degrees C/km)
df$spatg <- (df$ftas - df$tas)/df$dist

# Get the temporal gradient (degrees C/yr)
df$vector <- df$tempg/df$spatg



### PLOT VECTOR FIELD MAP
proj <- CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000
+ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")
vdf <- select(df, x, y, vector)

r <- rasterFromXYZ(vdf, crs=proj)
vectorplot(r, par.settings=RdBuTheme())
