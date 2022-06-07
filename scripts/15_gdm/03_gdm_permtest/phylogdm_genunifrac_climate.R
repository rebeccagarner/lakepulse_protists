rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_protists_genunifrac_climate.rda")
phylogdm_climate <- gdm.varImp(sitepair_climate_phylo, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 2)
save(phylogdm_climate, file = "lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_climate_100perm.rda")
