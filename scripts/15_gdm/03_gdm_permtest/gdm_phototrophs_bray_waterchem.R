rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_phototrophs_bray_waterchem.rda")
gdm_phototrophs_waterchem <- gdm.varImp(phototrophs_sitepair_waterchem, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 8)
save(gdm_phototrophs_waterchem, file = "lp2017-2019watercolumn18s_gdm_phototrophs_bray_waterchem_100perm.rda")
