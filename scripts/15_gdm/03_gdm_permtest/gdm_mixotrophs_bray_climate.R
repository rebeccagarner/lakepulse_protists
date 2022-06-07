rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_climate.rda")
gdm_mixotrophs_climate <- gdm.varImp(mixotrophs_sitepair_climate, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 4)
save(gdm_mixotrophs_climate, file = "lp2017-2019watercolumn18s_gdm_mixotrophs_bray_climate_100perm.rda")
