rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_phototrophs_bray_climate.rda")
gdm_phototrophs_climate <- gdm.varImp(phototrophs_sitepair_climate, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 4)
save(gdm_phototrophs_climate, file = "lp2017-2019watercolumn18s_gdm_phototrophs_bray_climate_100perm.rda")
