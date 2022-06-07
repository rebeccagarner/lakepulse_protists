rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_consumers_bray_climate.rda")
gdm_consumers_climate <- gdm.varImp(consumers_sitepair_climate, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 4)
save(gdm_consumers_climate, file = "lp2017-2019watercolumn18s_gdm_consumers_bray_climate_100perm.rda")
