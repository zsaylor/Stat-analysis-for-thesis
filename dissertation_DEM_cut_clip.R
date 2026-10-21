# incarcare EU DEM


library(terra)
library(sf)

############################################################
## DEM → EPSG:3844 → CLIP ROMANIA
############################################################

library(terra)
library(sf)

############################
## 1. LOAD DATA
############################

# DEM in EPSG:4326
dem <- rast("C:/Users/user/Documents/10_DEM_y40x20/10_DEM_y40x20.tif")

# Romania border (shp)
romania <- vect("F:/01_Zane/ro_admin_country_polygon/ro_admin_country_polygon.shp")


############################
## 2. REPROJECT DEM → EPSG:3844
##
## Important:
## - bilinear for continuous data (DEM)
############################

dem_3844 <- project(
  dem,
  "EPSG:3844",
  method = "bilinear"
)

############################
## 3. ENSURE ROMANIA CRS MATCHES DEM
############################

romania <- project(romania, crs(dem_3844))

############################
## 4. CROP + MASK (CLIP)
############################

dem_crop <- crop(dem_3844, romania)
dem_clip <- mask(dem_crop, romania)

############################
## 5. SAVE RESULT
############################

writeRaster(
  dem_clip,
  "DEM_Romania_EPSG3844.tif",
  overwrite = TRUE
)

############################
## 6. QUICK CHECK
############################

plot(dem_clip, col = terrain.colors(100), main = "DEM Romania (EPSG:3844)")

############################################################
## DEM CLASSIFICATION INTO ALTITUDE CLASSES
############################################################

# Reclassification matrix: from, to, new class
rcl <- matrix(c(
  -Inf,  200, 1,
  200,  500, 2,
  500,  800, 3,
  800, 1500, 4,
  1500,  Inf, 5
), ncol = 3, byrow = TRUE)

dem_class <- classify(dem_clip, rcl)


dem_class <- as.factor(dem_class)

levels(dem_class) <- data.frame(
  value = 1:5,
  label = c(
    "Plains (<200 m)",
    "Tablelands / Hills (200–500 m)",
    "Low Mountains (500–800 m)",
    "Medium-High Mountains (800–1500 m)",
    "Alpine (>1500 m)"
  )
)

cols <- c(
  "#4daf4a",  # plains – green
  "#a1d76a",  # hills – light green
  "#ffffbf",  # high hills – yellow
  "#fdae61",  # mountains – orange
  "#d73027"   # alpine – red
)

plot(
  dem_class,
  col = cols,
  main = "Elevation Zones – Romania"
)