rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_protists_genunifrac_morphometry.rda")
phylogdm_morphometry <- gdm.varImp(sitepair_morphometry_phylo, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 6)
save(phylogdm_morphometry, file = "lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_morphometry_100perm.rda")
