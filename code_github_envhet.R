#### Code accompanying: Schmidt, Dray, & Garroway 
# Genetic and species-level biodiversity patterns are linked by demography and ecological opportunity
# https://onlinelibrary.wiley.com/doi/abs/10.1111/evo.14407
# Dryad data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.g79cnp5qv

library(tidyverse)
library(sf)
library(raster)

# Code to extract environmental heterogeneity within 4 area sizes around sites (5k, 20k, 50k, 100k km2)
# Note range-wide mean PET and AET were extracted using species distribution shapefiles from the IUCN in ArcMap

## Data: 
sgdata <- read.csv("Schmidt_genetic_environments_dryad.csv", head = T)

# Convert to sf object:
sites <- st_as_sf(sgdata, coords = c("lon", "lat"), crs=4326)

# Landcover dataset - CEC (250m): 
# http://www.cec.org/north-american-environmental-atlas/land-cover-2010-modis-250m/
landcov250 <- raster("Land_Cover_2010_TIFF/LandCover_2010/data/NA_LandCover_2010_25haMMU.tif")

# Reproject sites into same coordinate system as raster:
sites_laea <- st_transform(sites, crs = crs(landcov))

## Simpson's index
simpson_diversity <- function(cat.vect){
  px  = table(cat.vect)/length(cat.vect)
  D1 = 1-sum(px^2)
  return(D1)
}

## 5k buffer
lc_5ke <- raster::extract(landcov, sites_laea, buffer = 40*1000) # outputs a list, every item is the summary for 1 site
lc_5kSI <- lapply(lc_5ke, function(x) simpson_diversity(as.factor(x)))
sgdata$lc_5kSI <- unlist(lc_5kSI)
rm(lc_5ke, lc_5kSI)

## 20k buffer
lc_20ke <- raster::extract(landcov, sites_laea, buffer = 80*1000)
lc_20kSI <- lapply(lc_20ke, function(x) simpson_diversity(as.factor(x)))
sgdata$lc_20kSI <- unlist(lc_20kSI)
rm(lc_20ke)

## 50k buffer
lc_50ke <- raster::extract(landcov, sites_laea, buffer = 126*1000)
lc_50k <- lapply(lc_50ke, function(x) length(unique(x)))
sgdata$lc_50k <- unlist(lc_50k)
rm(lc_50ke)

## 100k buffer
lc_100ke <- raster::extract(landcov, sites_laea, buffer = 178*1000)
lc_100k <- lapply(lc_100ke, function(x) length(unique(x)))
sgdata$lc_100k <- unlist(lc_100k)
rm(lc_100ke)