# Compute spatial medians within the first two dimensions of the ordination space

# Load libraries
library(tidyverse)
library(vegan)
library(ape)
library(janitor)

# Load palettes
source("scripts/00_palettes.R")


#### Import data ####
# Import (generalized alpha = 0.5) UniFrac distance matrix
load("output/dissim/lp2017-2019watercolumn18s_protists_unifrac_gen0pt5alpha.rda")

# Import melted protist dataset
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_function.tsv", col_names = TRUE)

# Import metadata
source("scripts/05_metadata.R")


#### Format ASV tables ####
# Create ASV table
asvtable <- asv_melt %>%
  select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id")
dim(asvtable)

# Divide ASVs by trophic functional group
# Recalculate relative sequence abundances
# Phototrophs ----
phototrophs_melt <- asv_melt %>%
  filter(trophic_group == "phototrophs") %>%
  dplyr::select(-relseqs)

phototrophs_nseqs <- phototrophs_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

phototrophs_melt <- phototrophs_melt %>%
  left_join(phototrophs_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  dplyr::select(-nseqs_total)
sum(phototrophs_melt$relseqs)  # Should evaluate to the number of lakes

asvtable_phototrophs <- phototrophs_melt %>%
  select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  column_to_rownames("lake_id")

# Consumers ----
consumers_melt <- asv_melt %>%
  filter(trophic_group == "consumers") %>%
  dplyr::select(-relseqs)

consumers_nseqs <- consumers_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

consumers_melt <- consumers_melt %>%
  left_join(consumers_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  dplyr::select(-nseqs_total)
sum(consumers_melt$relseqs)  # Should evaluate to the number of lakes

asvtable_consumers <- consumers_melt %>%
  select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  column_to_rownames("lake_id")

# Mixotrophs ----
mixotrophs_melt <- asv_melt %>%
  filter(trophic_group == "mixotrophs") %>%
  dplyr::select(-relseqs)

mixotrophs_nseqs <- mixotrophs_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

mixotrophs_melt <- mixotrophs_melt %>%
  left_join(mixotrophs_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  dplyr::select(-nseqs_total)
sum(mixotrophs_melt$relseqs)  # Should evaluate to the number of lakes

asvtable_mixotrophs <- mixotrophs_melt %>%
  select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  column_to_rownames("lake_id")


#### Compute spatial medians ####
# Compute distances
asv_bray <- vegdist(asvtable, method = "bray")
phototrophs_bray <- vegdist(asvtable_phototrophs, method = "bray")
consumers_bray <- vegdist(asvtable_consumers, method = "bray")
mixotrophs_bray <- vegdist(asvtable_mixotrophs, method = "bray")

# Define groups (trophic state)
metadata <- metadata %>%
  arrange(lake_id)
groups_trophic <- factor(metadata$trophic_status)

trophic <- metadata %>%
  dplyr::select(lake_id, trophic_status)

# Compute dispersions
tax_dispersions <- betadisper(asv_bray, group = groups_trophic, type = "median", sqrt.dist = FALSE)
phylo_dispersions <- betadisper(unifrac_gen0.5alpha_dist, group = groups_trophic, type = "median", sqrt.dist = FALSE)

phototrophs_dispersions <- betadisper(phototrophs_bray, group = groups_trophic, type = "median", sqrt.dist = FALSE)
consumers_dispersions <- betadisper(consumers_bray, group = groups_trophic, type = "median", sqrt.dist = FALSE)
mixotrophs_dispersions <- betadisper(mixotrophs_bray, group = groups_trophic, type = "median", sqrt.dist = FALSE)

# Define function to compute spatial medians of the first two PCoA coordinates
computeDistMedianPCoA1to2 <- function(dispersions) {
  pcoa_sites <- dispersions$vectors %>%
    as_tibble(rownames = "lake_id") %>%
    select(lake_id, PCoA1, PCoA2)
  
  centroids_pcoa1to2 <- dispersions$centroids %>%
    as_tibble(rownames = "trophic_status") %>%
    select(trophic_status, PCoA1, PCoA2) %>%
    rename(centroid_PCoA1 = PCoA1,
           centroid_PCoA2 = PCoA2)
  
  dist_median_pcoa1to2 <- pcoa_sites %>%
    left_join(trophic, by = "lake_id") %>%
    left_join(centroids_pcoa1to2, by = "trophic_status") %>%
    mutate(dist_median_pcoa1to2 = sqrt((PCoA1 - centroid_PCoA1)^2 + (PCoA2 - centroid_PCoA2)^2))
  
  return(dist_median_pcoa1to2)
}

tax_dist_median_pcoa1to2 <- computeDistMedianPCoA1to2(tax_dispersions) %>%
  mutate(dist_centroid_metric = "taxonomic")
phylo_dist_median_pcoa1to2 <- computeDistMedianPCoA1to2(phylo_dispersions) %>%
  mutate(dist_centroid_metric = "phylogenetic")

phototrophs_dist_median_pcoa1to2 <- computeDistMedianPCoA1to2(phototrophs_dispersions) %>%
  mutate(dist_centroid_metric = "phototrophs")
consumers_dist_median_pcoa1to2 <- computeDistMedianPCoA1to2(consumers_dispersions) %>%
  mutate(dist_centroid_metric = "consumers")
mixotrophs_dist_median_pcoa1to2 <- computeDistMedianPCoA1to2(mixotrophs_dispersions) %>%
  mutate(dist_centroid_metric = "mixotrophs")

# Plot box plots of distances to spatial medians
(dist_median_boxplots <- bind_rows(tax_dist_median_pcoa1to2,
                                   phylo_dist_median_pcoa1to2) %>%
    left_join(ecozones, by = "lake_id") %>%
    ggplot() +
    facet_wrap(~factor(dist_centroid_metric, levels = c("taxonomic", "phylogenetic")),
               ncol = 1) +
    geom_boxplot(aes(x = factor(trophic_status, levels = c("ultraoligotrophic", "oligotrophic", "mesotrophic", "mesoeutrophic", "eutrophic", "hypereutrophic")),
                     y = dist_median_pcoa1to2,
                     colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
                 outlier.shape = NA) +
    geom_jitter(aes(x = factor(trophic_status, levels = c("ultraoligotrophic", "oligotrophic", "mesotrophic", "mesoeutrophic", "eutrophic", "hypereutrophic")),
                    y = dist_median_pcoa1to2,
                    shape = factor(ecozone, levels = ecozone_longitude_order),
                    colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
                alpha = 0.7, size = 2, height = 0) +
    labs(y = "Distance from spatial median") +
    coord_flip() +
    scale_colour_manual(values = palette_trophic, name = "Trophic state") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    theme_light() %+replace%
    theme(strip.background = element_rect(fill = "black"), strip.text = element_text(colour = "white", face = "bold"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()))

(dist_median_functions_boxplots <- bind_rows(phototrophs_dist_median_pcoa1to2,
                                             consumers_dist_median_pcoa1to2,
                                             mixotrophs_dist_median_pcoa1to2) %>%
    left_join(ecozones, by = "lake_id") %>%
    ggplot() +
    facet_wrap(~factor(dist_centroid_metric, levels = c("phototrophs", "consumers", "mixotrophs")),
               ncol = 1) +
    geom_boxplot(aes(x = factor(trophic_status, levels = c("ultraoligotrophic", "oligotrophic", "mesotrophic", "mesoeutrophic", "eutrophic", "hypereutrophic")),
                     y = dist_median_pcoa1to2,
                     colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
                 outlier.shape = NA) +
    geom_jitter(aes(x = factor(trophic_status, levels = c("ultraoligotrophic", "oligotrophic", "mesotrophic", "mesoeutrophic", "eutrophic", "hypereutrophic")),
                    y = dist_median_pcoa1to2,
                    shape = factor(ecozone, levels = ecozone_longitude_order),
                    colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
                alpha = 0.7, size = 2, height = 0) +
    labs(y = "Distance from spatial median") +
    coord_flip() +
    scale_colour_manual(values = palette_trophic, name = "Trophic state") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    theme_light() %+replace%
    theme(strip.background = element_rect(fill = "black"), strip.text = element_text(colour = "white", face = "bold"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()))
#ggsave("figures/lp2017-2019watercolumn18s_trophic_dist_median_functions.pdf", dist_median_functions_boxplots, "pdf", width = 8, height = 9, units = "in")


#### Compare distances to spatial medians ####
# Define function to perform Tukey's test
testTukeyDistMedian <- function(dispersions) {
  pcoa_sites <- dispersions$vectors %>%
    as_tibble(rownames = "lake_id") %>%
    select(lake_id, PCoA1, PCoA2)
  
  centroids_pcoa1to2 <- dispersions$centroids %>%
    as_tibble(rownames = "trophic_status") %>%
    select(trophic_status, PCoA1, PCoA2) %>%
    rename(centroid_PCoA1 = PCoA1,
           centroid_PCoA2 = PCoA2)
  
  dist_median_pcoa1to2 <- pcoa_sites %>%
    left_join(trophic, by = "lake_id") %>%
    left_join(centroids_pcoa1to2, by = "trophic_status") %>%
    mutate(dist_median_pcoa1to2 = sqrt((PCoA1 - centroid_PCoA1)^2 + (PCoA2 - centroid_PCoA2)^2))
  
  anova_dist_median_pcoa1to2 <- aov(dist_median_pcoa1to2 ~ trophic_status, data = dist_median_pcoa1to2)
  
  tukey_dist_median_pcoa1to2 <- TukeyHSD(anova_dist_median_pcoa1to2)$trophic_status %>%
    as_tibble(rownames = "comparison") %>%
    clean_names() %>%
    filter(!grepl("ultraoligotrophic", comparison)) %>%
    #filter(p_adj < 0.05) %>%
    mutate(significance_level = case_when(p_adj < 0.05 & p_adj >= 0.01 ~ "*",
                                          p_adj < 0.01 & p_adj >= 0.001 ~ "**",
                                          p_adj < 0.001 ~ "***",
                                          TRUE ~ "NS")) %>%
    rename(difference_in_means = diff,
           lower_endpoint_interval = lwr,
           upper_endpoint_interval = upr)
  
  return(tukey_dist_median_pcoa1to2)
}

tax_dist_median_tukey <- testTukeyDistMedian(tax_dispersions) %>%
  mutate(dist_median_metric = "taxonomic")
phylo_dist_median_tukey <- testTukeyDistMedian(phylo_dispersions) %>%
  mutate(dist_median_metric = "phylogenetic")

phototrophs_dist_median_tukey <- testTukeyDistMedian(phototrophs_dispersions) %>%
  mutate(dist_median_metric = "phototrophs")
consumers_dist_median_tukey <- testTukeyDistMedian(consumers_dispersions) %>%
  mutate(dist_median_metric = "consumers")
mixotrophs_dist_median_tukey <- testTukeyDistMedian(mixotrophs_dispersions) %>%
  mutate(dist_median_metric = "mixotrophs")

dist_median_tukey_all <- bind_rows(tax_dist_median_tukey,
          phylo_dist_median_tukey,
          phototrophs_dist_median_tukey,
          consumers_dist_median_tukey,
          mixotrophs_dist_median_tukey) %>%
  filter(p_adj < 0.05) %>%
  select(dist_median_metric, comparison, difference_in_means, lower_endpoint_interval, upper_endpoint_interval, p_adj, significance_level)

# dist_median_tukey_all %>%
#   write_csv("output/diversity/lp2017-2019watercolumn18s_trophic_dist_median_tukey.csv")
