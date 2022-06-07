rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_consumers_bray_landuse_soil.rda")
gdm_consumers_landuse_soil <- gdm.varImp(consumers_sitepair_landuse_soil, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 10)
save(gdm_consumers_landuse_soil, file = "lp2017-2019watercolumn18s_gdm_consumers_bray_landuse_soil_100perm.rda")
