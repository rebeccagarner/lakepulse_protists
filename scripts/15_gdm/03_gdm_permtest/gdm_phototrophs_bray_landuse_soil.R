rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_phototrophs_bray_landuse_soil.rda")
gdm_phototrophs_landuse_soil <- gdm.varImp(phototrophs_sitepair_landuse_soil, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 10)
save(gdm_phototrophs_landuse_soil, file = "lp2017-2019watercolumn18s_gdm_phototrophs_bray_landuse_soil_100perm.rda")
