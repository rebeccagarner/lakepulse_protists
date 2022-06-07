# Compute ASV rarefaction across the landscape (all lakes)

# Load libraries
library(tidyverse)
library(vegan)


#### Import and format data ####
# Import melted protists data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)


#### Create taxonomy table ####
taxonomy <- asv_melt %>%
  distinct(asv_code, kingdom, supergroup, division, class, order, family, genus, species)


#### Collapse all ASVs into one landscape-wide pooled ASV table ####
asvtable_landscape <- asv_melt %>%
  group_by(asv_code) %>%
  summarize(nseqs = sum(nseqs)) %>%
  pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0) %>%
  as.data.frame()
rownames(asvtable_landscape) <- c("landscape")


#### Rarefaction curve ####
rarefaction <- tibble(rarefied = 0, richness = 0)
rarefied <- 1
while (rarefied <= sum(asvtable_landscape)) {
  asvtable_rarefied <- rrarefy(asvtable_landscape, sample = rarefied)
  richness <- specnumber(asvtable_rarefied)
  rarefaction_tmp <- tibble(rarefied, richness)
  rarefaction <- bind_rows(rarefaction, rarefaction_tmp)
  rarefied <- rarefied + 1000
  print(rarefied)
}

rarefaction_allseqs <- tibble(rarefied = sum(asvtable_landscape), richness = specnumber(asvtable_landscape))
rarefaction <- bind_rows(rarefaction, rarefaction_allseqs)

# rarefaction %>%
#   write_tsv("output/diversity/lp2017-2019watercolumn18s_protists_rarefaction.tsv", col_names = TRUE)

(rarefaction_curve <- rarefaction %>%
    ggplot(aes(x = rarefied, y = richness)) +
    geom_line() +
    labs(x = "Sample size (sequences)", y = "ASV richness") +
    theme_light())
#ggsave(filename = "figures/lp2017-2019watercolumn18s_protists_rarefaction_curve.pdf", plot = rarefaction_curve, device = "pdf", width = 6, height = 5, units = "in")
