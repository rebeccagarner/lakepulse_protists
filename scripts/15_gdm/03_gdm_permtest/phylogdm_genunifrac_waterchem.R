rm(list = ls())
library(gdm)

load("lp2017-2019watercolumn18s_sitepair_protists_genunifrac_waterchem.rda")
phylogdm_waterchem <- gdm.varImp(sitepair_waterchem_phylo, geo = FALSE, nPerm = 100, fullModelOnly = FALSE, parallel = TRUE, cores = 8)
save(phylogdm_waterchem, file = "lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_waterchem_100perm.rda")
