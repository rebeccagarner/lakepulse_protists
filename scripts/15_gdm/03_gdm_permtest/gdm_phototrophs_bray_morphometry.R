rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_phototrophs_bray_morphometry.rda")
gdm_phototrophs_morphometry <- gdm.varImp(phototrophs_sitepair_morphometry, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 6)
save(gdm_phototrophs_morphometry, file = "lp2017-2019watercolumn18s_gdm_phototrophs_bray_morphometry_100perm.rda")
