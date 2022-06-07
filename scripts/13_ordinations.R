# Perform ordinations (PCA, PCoA)

# Load libraries
library(tidyverse)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ade4)
library(ape)
library(ggpubr)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import melted protist dataset
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_function.tsv", col_names = TRUE)

# Import metadata
source("scripts/05_metadata.R")

# Import phylogenetic (UniFrac) distance matrices
load("output/dissim/lp2017-2019watercolumn18s_protists_unifrac_gen0pt5alpha.rda")


#### Create taxonomy table ####
taxonomy <- asv_melt %>%
  distinct(asv_code, kingdom, supergroup, division, class, order, family, genus, species)


#### Create community data table ####
# ASV table based on protist dataset
asvtable <- asv_melt %>%
  select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  column_to_rownames("lake_id")

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

# Parasites ----
parasites_melt <- asv_melt %>%
  filter(trophic_group == "parasites" | trophic_group == "consumers, parasites") %>%
  dplyr::select(-relseqs)

parasites_nseqs <- parasites_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

parasites_melt <- parasites_melt %>%
  left_join(parasites_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  dplyr::select(-nseqs_total)
sum(parasites_melt$relseqs)  # Should evaluate to the number of lakes

asvtable_parasites <- parasites_melt %>%
  select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  column_to_rownames("lake_id")


#### Transform community data ####
# Hellinger transformation
asvtable_hellinger <- decostand(asvtable, "hellinger")

asvtable_phototrophs_hellinger <- decostand(asvtable_phototrophs, "hellinger")
asvtable_consumers_hellinger <- decostand(asvtable_consumers, "hellinger")
asvtable_mixotrophs_hellinger <- decostand(asvtable_mixotrophs, "hellinger")
asvtable_parasites_hellinger <- decostand(asvtable_parasites, "hellinger")


#### PCA ####
# Define function to plot PCA with dark theme
plotPCA <- function(community_data, arrows_taxa = 1) {
  # Compute PCA
  pca <- rda(community_data)
  summary(pca)  # Summarize PCA results
  
  # Extract PCA site and species scores
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "lake_id") %>%
    left_join(metadata, by = "lake_id")
  pca_species <- scores(pca)$species %>%
    as_tibble(rownames = "asv_code")
  
  # Calculate variation explained by PCA dimensions 1 and 2
  pc1_pctvar <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1)
  pc2_pctvar <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1)
  
  # PC contributions
  contributors <- pca_species %>%
    mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
    slice_max(magnitude, n = 7) %>%
    pull(asv_code)
  
  # Curate taxon contributor labels in PCA plot
  contribution <- pca_species %>%
    filter(asv_code %in% contributors) %>%
    left_join(taxonomy, by = "asv_code") %>%
    mutate(label = case_when(!is.na(species) ~ str_c(species, "\n(", division, ")", sep = ""),
                             is.na(species) & !is.na(genus) ~ str_c(genus, "\n(", division, ")", sep = ""),
                             is.na(genus) & !is.na(family) ~ str_c(family, "\n(", division, ")", sep = ""),
                             is.na(family) & !is.na(order) ~ str_c(order, "\n(", division, ")", sep = ""),
                             is.na(order) & !is.na(class) ~ str_c(class, "\n(", division, ")", sep = ""),
                             is.na(class) & !is.na(division) ~ str_c(division),
                             is.na(division) & !is.na(supergroup) ~ str_c(supergroup))) %>%
    mutate(label = str_replace_all(label, "_", " "))
  
  # Plot PCA
  pca_plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#E6E6E6") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#E6E6E6") +
    geom_text(aes(x = (arrows_taxa + 0.1)*PC1, y = (arrows_taxa + 0.1)*PC2, label = label),
              data = contribution,
              size = 5, colour = "#F8F8F8", fontface = "italic") +
    geom_segment(aes(x = 0, y = 0, xend = arrows_taxa*PC1, yend = arrows_taxa*PC2),
                 data = contribution,
                 arrow = arrow(length = unit(x = 0.2, units = "cm")),
                 alpha = 0.7, colour = "#A0A0A0") +
    geom_point(aes(x = PC1, y = PC2,
                   colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic")),
                   shape = factor(ecozone, levels = ecozone_longitude_order)),
               data = pca_sites,
               size = 5, alpha = 0.9, stroke = 2) +
  scale_colour_manual(values = palette_trophic) +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    labs(x = paste("PC axis 1 (", pc1_pctvar, "%)", sep = ""),
         y = paste("PC axis 2 (", pc2_pctvar, "%)", sep = ""),
         colour = "Trophic status") +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "#383838"),
          legend.key = element_blank(), legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  return(pca_plot)
}

# PCA of Hellinger-transformed ASV data
(pca_hellinger <- plotPCA(asvtable_hellinger, arrows_taxa = 1.8))
#ggsave("figures/lp2017-2019watercolumn18s_hellinger_pca_dark.pdf", pca_hellinger, "pdf", width = 12, height = 12, units = "in")

# PCA on Hellinger-transformed community data by trophic mode
pca_phototrophs_hellinger <- plotPCA(asvtable_phototrophs_hellinger, arrows_taxa = 1.4)
pca_consumers_hellinger <- plotPCA(asvtable_consumers_hellinger, arrows_taxa = 0.8)
pca_mixotrophs_hellinger <- plotPCA(asvtable_mixotrophs_hellinger, arrows_taxa = 0.4)
(pca_hellinger_all <- ggarrange(pca_phototrophs_hellinger, pca_consumers_hellinger, pca_mixotrophs_hellinger,
                                ncol = 1, common.legend = TRUE, legend = "none"))
#ggsave("figures/lp2017-2019watercolumn18s_functions_hellinger_pca_dark.pdf", pca_hellinger_all, "pdf", width = 12, height = 36, units = "in")


#### PCoA ####
plotPCoA <- function(dist) {
  # Compute PCoA
  prcoorda <- pcoa(dist)
  
  # Extract site coordinates
  pcoa_sites <- prcoorda$vectors[,c("Axis.1", "Axis.2")] %>%
    as_tibble(rownames = "lake_id") %>%
    left_join(metadata, by = "lake_id")
  
  # Calculate variation explained by PCoA dimensions 1 and 2
  (pcoa1_pctvar <- round(prcoorda$values$Relative_eig[1] * 100, digits = 1))
  (pcoa2_pctvar <- round(prcoorda$values$Relative_eig[2] * 100, digits = 1))
  
  # Plot PCoA
  pcoa_plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#E6E6E6") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#E6E6E6") +
    geom_point(aes(x = Axis.1, y = Axis.2,
                   colour = factor(trophic_status,
                                   levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic")),
                   shape = factor(ecozone, levels = ecozone_longitude_order)),
               data = pcoa_sites,
               size = 5, alpha = 0.9, stroke = 2) +
    scale_colour_manual(values = palette_trophic, name = "Trophic status") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    labs(x = paste("PCoA axis 1 (", pcoa1_pctvar, "%)", sep = ""),
         y = paste("PCoA axis 2 (", pcoa2_pctvar, "%)", sep = "")) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "#383838"),
          legend.key = element_blank(), legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  return(pcoa_plot)
}

# PCoA on a generalized (alpha = 0.5) UniFrac distance matrix
(pcoa_unifrac_gen0.5alpha <- plotPCoA(unifrac_gen0.5alpha_dist))
#ggsave("figures/lp2017-2019watercolumn18s_genunifrac_pcoa_dark.pdf", pcoa_unifrac_gen0.5alpha, "pdf", width = 12, height = 12, units = "in")


#### Compute RV coefficient ####
# Compute PCA
asvs_hellinger_pca <- rda(asvtable_hellinger)
summary(asvs_hellinger_pca)

# Extract PCA site coordinates
asvs_hellinger_pca_sites1to2 <- scores(asvs_hellinger_pca, choices = c(1,2))$sites %>%
  as_tibble(rownames = "lake_id") %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id")

# Compute PCoA
genunifrac_pcoa <- pcoa(unifrac_gen0.5alpha_dist)

# Extract PCoA site coordinates
genunifrac_pcoa_sites1to2 <- genunifrac_pcoa$vectors[,c("Axis.1", "Axis.2")] %>%
  as_tibble(rownames = "lake_id") %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id")

# Compute RV coefficient
coeffRV(asvs_hellinger_pca_sites1to2, genunifrac_pcoa_sites1to2)
