rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_protists_genunifrac_geography.rda")
phylogdm_geography <- gdm.varImp(sitepair_geography_phylo, geo = TRUE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 2)
save(phylogdm_geography, file = "lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_geography_100perm.rda")
