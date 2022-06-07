rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_consumers_bray_waterchem.rda")
gdm_consumers_waterchem <- gdm.varImp(consumers_sitepair_waterchem, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 8)
save(gdm_consumers_waterchem, file = "lp2017-2019watercolumn18s_gdm_consumers_bray_waterchem_100perm.rda")
