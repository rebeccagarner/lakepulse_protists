rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_protists_genunifrac_landuse_soil.rda")
phylogdm_landuse_soil <- gdm.varImp(sitepair_landuse_soil_phylo, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 10)
save(phylogdm_landuse_soil, file = "lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_landuse_soil_100perm.rda")
