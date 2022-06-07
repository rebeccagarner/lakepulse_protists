# Create phylogenetic dissimilarity matrices

# Load libraries
library(tidyverse)
library(vegan)
library(ape)
library(GUniFrac)
library(picante)


#### Import and format data ####
# Import rooted phylogenetic tree
tree <- read.tree("output/trees/lp2017-2019watercolumn18s_protists_366lakes_sina_gtr_rooted.tree")

# Import melted protist ASV dataset
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv")
length(unique(asv_melt$asv_code)) == length(tree$tip.label)  # Should evaluate to TRUE

# Format ASV table
asvtable <- asv_melt %>%
  select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  column_to_rownames("lake_id")

# Match phylogeny and ASV table
asv_tree_matched <- match.phylo.comm(phy = tree, comm = asvtable)


#### Calculate UniFrac distances ####
# Calculate weighted, unweighted, and generalized UniFrac distances
asv_unifrac <- GUniFrac(asv_tree_matched$comm, asv_tree_matched$phy, alpha = c(0, 0.5, 1))$unifracs

unifrac_weighted <- asv_unifrac[, , "d_1"]  # Weighted UniFrac
unifrac_weighted_dist <- as.dist(unifrac_weighted)
#save(unifrac_weighted_dist, file = "output/dissim/lp2017-2019watercolumn18s_protists_unifrac_weighted.rda")
unifrac_unweighted <- asv_unifrac[, , "d_UW"]  # Unweighted UniFrac
unifrac_unweighted_dist <- as.dist(unifrac_unweighted)
#save(unifrac_unweighted_dist, file = "output/dissim/lp2017-2019watercolumn18s_protists_unifrac_unweighted.rda")
unifrac_gen0.5alpha <- asv_unifrac[, , "d_0.5"]  # Generalized UniFrac with alpha 0.5
unifrac_gen0.5alpha_dist <- as.dist(unifrac_gen0.5alpha)
#save(unifrac_gen0.5alpha_dist, file = "output/dissim/lp2017-2019watercolumn18s_protists_unifrac_gen0pt5alpha.rda")


#### Plot pair-wise UniFrac distances ####
# Define function to return upper triangle of distance matrix
halveDist <- function(mat) {
  for (i in 1:ncol(mat)) {
    for (j in 1:i) {
      mat[i,j] <- NA
    }
  }
  
  return(mat)
}

unifrac_weighted_unpaired <- halveDist(unifrac_weighted) %>%
  as_tibble(rownames = "lake_id1") %>%
  pivot_longer(!lake_id1, names_to = "lake_id2", values_to = "dist") %>%
  filter(!is.na(dist)) %>%
  mutate(unifrac = "weighted")

unifrac_unweighted_unpaired <- halveDist(unifrac_unweighted) %>%
  as_tibble(rownames = "lake_id1") %>%
  pivot_longer(!lake_id1, names_to = "lake_id2", values_to = "dist") %>%
  filter(!is.na(dist)) %>%
  mutate(unifrac = "unweighted")

unifrac_gen0.5alpha_unpaired <- halveDist(unifrac_gen0.5alpha) %>%
  as_tibble(rownames = "lake_id1") %>%
  pivot_longer(!lake_id1, names_to = "lake_id2", values_to = "dist") %>%
  filter(!is.na(dist)) %>%
  mutate(unifrac = "gen0.5alpha")

(unifrac_hist <- bind_rows(unifrac_weighted_unpaired,
                           unifrac_unweighted_unpaired,
                           unifrac_gen0.5alpha_unpaired) %>%
    ggplot() +
    facet_wrap(~factor(unifrac, levels = c("weighted", "gen0.5alpha", "unweighted"))) +
    geom_histogram(aes(dist), binwidth = 0.01) +
    labs(x = "UniFrac distance",
         y = "Number of site-pairs") +
    theme_bw())
#ggsave("figures/lp2017-2019watercolumn18s_unifrac_histograms.pdf", unifrac_hist, "pdf", width = 10, height = 5, units = "in")
