# Calculate geographic distances for the study

# Load libraries
library(tidyverse)
library(GeoRange)


#### Import and format data ####
# Import metadata
source("scripts/05_metadata.R")


#### Calculate geographic ranges ####
# Calculate convex hull area in km^2 (i.e. polygon area enclosed by outermost points) from lake coordinates
CHullAreaEarth(longs = metadata$longitude, lats = metadata$latitude)

# Calculate the latitudinal range in degrees and km
LatRg(metadata$latitude)

# Calculate the longitudinal range in degrees and km
LongRg(metadata$longitude)
