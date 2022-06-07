rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_consumers_bray_morphometry.rda")
gdm_consumers_morphometry <- gdm.varImp(consumers_sitepair_morphometry, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 6)
save(gdm_consumers_morphometry, file = "lp2017-2019watercolumn18s_gdm_consumers_bray_morphometry_100perm.rda")
