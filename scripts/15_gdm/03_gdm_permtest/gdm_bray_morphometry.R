rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_protists_bray_morphometry.rda")
gdm_morphometry <- gdm.varImp(sitepair_morphometry, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 8)
save(gdm_morphometry, file = "lp2017-2019watercolumn18s_gdm_protists_bray_morphometry_100perm.rda")
