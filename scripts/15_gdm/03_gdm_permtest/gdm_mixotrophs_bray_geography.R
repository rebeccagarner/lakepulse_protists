rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_geography.rda")
gdm_mixotrophs_geography <- gdm.varImp(mixotrophs_sitepair_geography, geo = TRUE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 2)
save(gdm_mixotrophs_geography, file = "lp2017-2019watercolumn18s_gdm_mixotrophs_bray_geography_100perm.rda")
