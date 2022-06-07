rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_phototrophs_bray_geography.rda")
gdm_phototrophs_geography <- gdm.varImp(phototrophs_sitepair_geography, geo = TRUE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 2)
save(gdm_phototrophs_geography, file = "lp2017-2019watercolumn18s_gdm_phototrophs_bray_geography_100perm.rda")
