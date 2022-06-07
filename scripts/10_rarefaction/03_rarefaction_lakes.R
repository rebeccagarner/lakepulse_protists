# Compute ASV rarefaction curves by lake

# Load libraries
library(tidyverse)
library(vegan)


#### Import and format data ####
# Import melted protists data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)


#### Create taxonomy table ####
taxonomy <- asv_melt %>%
  distinct(asv_code, kingdom, supergroup, division, class, order, family, genus, species)


#### Create ASV tables ####
# ASV table of sequence counts
asvtable_nseqs <- asv_melt %>%
  select(lake_id, asv_code, nseqs) %>%
  spread(asv_code, nseqs, fill = 0) %>%
  column_to_rownames("lake_id")


#### Rarefaction curve ####
# Compute rarefaction curve for each site
for (i in 1:nrow(asvtable_nseqs)) {
  asvtable_row <- asvtable_nseqs[i,]
  lake_id <- rownames(asvtable_row)
  
  rarefaction <- tibble(rarefied = 0, richness = 0)
  rarefied <- 1
  while (rarefied <= rowSums(asvtable_row)) {
    asvtable_rarefied <- rrarefy(asvtable_row, sample = rarefied)
    richness <- specnumber(asvtable_rarefied)
    rarefaction_tmp <- tibble(rarefied, richness)
    rarefaction <- bind_rows(rarefaction, rarefaction_tmp)
    rarefied <- rarefied + 1000
    print(rarefied)
  }
  rarefaction_allseqs <- tibble(rarefied = sum(asvtable_row), richness = specnumber(asvtable_row))
  rarefaction <- bind_rows(rarefaction, rarefaction_allseqs)
  
  rarefaction %>%
    write_tsv(paste0("output/diversity/rarefaction/lp2017-2019watercolumn18s_protists_", lake_id, "_rarefaction.tsv"), col_names = TRUE)
}
