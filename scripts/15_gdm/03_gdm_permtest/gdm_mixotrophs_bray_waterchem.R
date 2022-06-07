rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_waterchem.rda")
gdm_mixotrophs_waterchem <- gdm.varImp(mixotrophs_sitepair_waterchem, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 8)
save(gdm_mixotrophs_waterchem, file = "lp2017-2019watercolumn18s_gdm_mixotrophs_bray_waterchem_100perm.rda")
