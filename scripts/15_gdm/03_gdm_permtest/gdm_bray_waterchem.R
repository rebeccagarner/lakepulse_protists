rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_protists_bray_waterchem.rda")
gdm_waterchem <- gdm.varImp(sitepair_waterchem, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 8)
save(gdm_waterchem, file = "lp2017-2019watercolumn18s_gdm_protists_bray_waterchem_100perm.rda")
