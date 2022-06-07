rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_protists_bray_climate.rda")
gdm_climate <- gdm.varImp(sitepair_climate, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 4)
save(gdm_climate, file = "lp2017-2019watercolumn18s_gdm_protists_bray_climate_100perm.rda")
