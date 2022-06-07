rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_morphometry.rda")
gdm_mixotrophs_morphometry <- gdm.varImp(mixotrophs_sitepair_morphometry, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 6)
save(gdm_mixotrophs_morphometry, file = "lp2017-2019watercolumn18s_gdm_mixotrophs_bray_morphometry_100perm.rda")
