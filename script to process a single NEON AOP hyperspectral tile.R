## ---------------------------
##
## Script name: Process a single NEON AOP hyperspectral file
##
## Purpose of script: This script will process a single NEON hyperspectral file, 
##                    similar to tutorials on their website. However this script
##                    will create an NDVI and shade mask that can be applied 
##                    to the processed tile 
##
## Author: Alexander Cox
##
## Email: acox22@nd.edu
## ----------------------------

#clear environment 
rm(list=ls())


# ------------------------------------------------------------------------------
# 1. Install and load required packages 
# ------------------------------------------------------------------------------


#install required packages 
list.of.packages <- c("BiocManager", "dplyr", "neonUtilities", "sf", "terra", 
                      "tidyr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#install rhdf5 from Bioconductor to work with the file type 
BiocManager::install("rhdf5")

#load required packages  
library(dplyr)
library(neonUtilities)
library(rhdf5)
library(sf)
library(terra)
library(tidyr)


# ------------------------------------------------------------------------------
# 2. Load in the DSM and corresponding hyperspectral file 
# ------------------------------------------------------------------------------


#load the digital surface model (DSM) file and set the coordinate reference 
#system. You may need to change it depending on location. This will be used to 
#create the shademask 
dsm <- rast("S:/Neon data/UNDE 2024/NEON_lidar-elev/NEON.D05.UNDE.DP3.30024.001.2024-07.basic.20250331T232419Z.PROVISIONAL/NEON_D05_UNDE_DP3_304000_5123000_DSM.tif")
crs(dsm) <- "EPSG:32616"

#load the hyperspectral file 
h5_file <- "F:/UNDE_2024-unprocessed-refl-tiles/NEON_refl-surf-bidir-ortho-mosaic/NEON.D05.UNDE.DP3.30006.002.2022-06.basic.20250114T202816Z.PROVISIONAL/NEON_D05_UNDE_DP3_304000_5123000_bidirectional_reflectance.h5"
wavelengthInfo <- h5readAttributes(h5_file,"/UNDE/Reflectance/Metadata/Spectral_Data/Wavelength")

#read in the wavelenghts 
wavelengths <- h5read(h5_file,"/UNDE/Reflectance/Metadata/Spectral_Data/Wavelength")

reflInfo <- h5readAttributes(h5_file, "/UNDE/Reflectance/Reflectance_Data")

nRows <- reflInfo$Dimensions[1]

nCols <- reflInfo$Dimensions[2]

nBands <- reflInfo$Dimensions[3]

h5EPSG <- h5read(h5_file, "/UNDE/Reflectance/Metadata/Coordinate_System/EPSG Code")

h5CRS <- crs(paste0("+init=epsg:",h5EPSG))

#CRS projection
new_crs_obj <- crs("+proj=utm +zone=16 +datum=WGS84")

#grab the UTM coordinates of the spatial extent
xMin <- reflInfo$Spatial_Extent_meters[1]
xMax <- reflInfo$Spatial_Extent_meters[2]
yMin <- reflInfo$Spatial_Extent_meters[3]
yMax <- reflInfo$Spatial_Extent_meters[4]

#define the extent (left, right, top, bottom)
rastExt <- ext(xMin,xMax,yMin,yMax)

##SpatExtent : 304000, 305000, 5123000, 5124000 (xmin, xmax, ymin, ymax)

#define the no data value
h5NoDataValue <- as.integer(reflInfo$Data_Ignore_Value)
cat('No Data Value:',h5NoDataValue)

#slice a band of data from the hdf5 file, extracts the reflectance array for that band,
#convert the data into a matrix, converts it to a raster, and then returns a spatially corrected raster
band2Raster <- function(file, band, noDataValue, extent, CRS){
  #first, read in the raster
  out <- h5read(file,"/UNDE/Reflectance/Reflectance_Data",index=list(band,NULL,NULL))
  #Convert from array to matrix
  out <- (out[1,,])
  #transpose data to fix flipped row and column order 
  #depending upon how your data are formatted you might not have to perform this
  #step.
  out <- t(out)
  #assign data ignore values to NA
  #note, you might chose to assign values of 15000 to NA
  out[out == noDataValue] <- NA
  
  #turn the out object into a raster
  outr <- rast(out,crs=as.character(CRS))
  
  #assign the extents to the raster
  ext(outr) <- extent
  
  #return the terra raster object
  return(outr)
}

#create a list of bands to include in the stack
rgb <- as.list(1:426)

rgb_rast <- lapply(rgb,FUN=band2Raster, file = h5_file,
                   noDataValue=h5NoDataValue, 
                   ext=rastExt,
                   CRS=h5CRS)

#create a raster stack from the list of rasters 
rgbStack <- rast(rgb_rast)

#create a list of band names
bandNames <- paste("Band_",unlist(rgb),sep="")

#create a vector of wavelength-based band names
bandNames <- paste("Wv", round(wavelengths, 2), sep = ".")

#set the rasterStack's names equal to the list of bandNames created above
names(rgbStack) <- bandNames

#scale the data as specified in the reflInfo$Scale Factor 
rgbStack <- rgbStack/as.integer(reflInfo$Scale_Factor)

head(rgbStack)

h5info <- h5ls(h5_file, recursive=TRUE)

#build a full path from group + name
h5info$fullpath <- file.path(h5info$group, h5info$name)

#locate azimuth/zenith by searching fullpath
az_idx <- grep("Solar_Azimuth_Angle", h5info$fullpath)
zn_idx <- grep("Solar_Zenith_Angle",  h5info$fullpath)

az_path <- h5info$fullpath[az_idx[1]]
zn_path <- h5info$fullpath[zn_idx[1]]

solar_az <- h5read(h5_file, az_path)
solar_zn <- h5read(h5_file, zn_path)

list(azimuth = solar_az, zenith = solar_zn)

makeShadeMask <- function(dsm, solar_az, solar_zn, threshold = 0.8) {
  #calculate slope and aspect in degrees 
  dsm_slope_deg  <- terrain(dsm, v = "slope",  unit="degrees")
  dsm_aspect_deg <- terrain(dsm, v = "aspect", unit="degrees")
  
  #calculate shade 
  s <- shade(
    slope     = dsm_slope_deg,
    aspect    = dsm_aspect_deg,
    angle     = solar_zn,  # degrees
    direction = solar_az,  # degrees
    normalize = TRUE
  )
  #set the threshold 
  s[s < threshold] <- NA
  return(s)
}

#apply the function, pixels with < 80% illumination will be masked
shade_mask <- makeShadeMask(dsm, solar_az, solar_zn, threshold = 0.8)

#make sure the resolution and extents match, resample if needed 
if (!compareGeom(rgbStack, shade_mask, stopOnError = FALSE)) {
  shade_mask <- resample(shade_mask, rgbStack, method="near")
}

#set shade masked pixels to NA 
rgbStack[is.na(shade_mask)] <- NA

#calculate NDVI from the bands in NEON's tutorial, apply a 0.8 threshold, 
#setting pixel values < 0.8 to NA 
red_band <- rgbStack[[58]]  # ~667 nm
nir_band <- rgbStack[[90]]  # ~827 nm
ndvi <- (nir_band - red_band) / (nir_band + red_band)
rgbStack[ndvi < 0.8] <- NA

#write out 
writeRaster(rgbStack, "S:/Neon data/test.masked.spectra.tif")