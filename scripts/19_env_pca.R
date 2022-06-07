# Perform ordinations of environmental data

# Load libraries
library(tidyverse)
library(vegan)
library(patchwork)
library(ggnewscale)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import metadata
source("scripts/05_metadata.R")


#### PCA ####
# Define function to plot PCA
plotPCA <- function(env_data, plot_title, arrows_taxa = 1) {
  # Format environmental data
  env_data <- env_data %>%
    filter(lake_id %in% lakes) %>%
    dplyr::select(lake_id, where(is.numeric)) %>%
    column_to_rownames("lake_id")
  
  # Scale and center data
  env_sc <- scale(env_data, center = TRUE, scale = TRUE)
  
  # Compute PCA
  pca <- rda(env_sc)
  #summary(pca)  # Summarize PCA results
  
  # Extract PCA site and species scores
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "lake_id") %>%
    left_join(metadata, by = "lake_id")
  pca_species <- scores(pca)$species %>%
    as_tibble(rownames = "variable") %>%
    left_join(variable_labels_nounits, by = "variable")
  
  # Calculate variation explained by PCA dimensions 1 and 2
  pc1_pctvar <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1)
  pc2_pctvar <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1)
  
  # Plot PCA
  pca_plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#E6E6E6") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#E6E6E6") +
    geom_text(aes(x = (arrows_taxa + 0.1)*PC1, y = (arrows_taxa + 0.1)*PC2, label = variable_curated),
              data = pca_species,
              colour = "#A0A0A0", size = 4) +
    geom_segment(aes(x = 0, y = 0, xend = arrows_taxa*PC1, yend = arrows_taxa*PC2),
                 data = pca_species,
                 arrow = arrow(length = unit(x = 0.2, units = "cm")),
                 alpha = 0.7, colour = "#A0A0A0") +
    geom_point(aes(x = PC1, y = PC2,
                   colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic")),
                   shape = factor(ecozone, levels = ecozone_longitude_order)),
               data = pca_sites, size = 5, alpha = 0.9, stroke = 2) +
    scale_colour_manual(values = palette_trophic) +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    labs(x = paste("PC axis 1 (", pc1_pctvar, "%)", sep = ""),
         y = paste("PC axis 2 (", pc2_pctvar, "%)", sep = ""),
         colour = "Trophic state") +
    ggtitle(plot_title) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          legend.key = element_blank(), legend.position = "right",
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  return(pca_plot)
}

# PCA on scaled and centered environmental data
(pca_geography <- plotPCA(geography %>%
                            dplyr::select(-contains("_x"), -contains("_y")), "Geography", arrows_taxa = 0.3))
(pca_morphometry <- plotPCA(morphometry, "Morphometry", arrows_taxa = 1.2))
(pca_landuse <- plotPCA(landuse, "Land use", arrows_taxa = 1.5))
(pca_soil <- plotPCA(soil, "Soil", arrows_taxa = 0.8))
(pca_climate <- plotPCA(climate, "Weather", arrows_taxa = 0.5))
(pca_waterchem <- plotPCA(waterchem, "Physicochemistry", arrows_taxa = 1.6))

(pca_env_all <- plotPCA(geography %>%
                          dplyr::select(-contains("_x"), -contains("_y"), -longitude) %>%
                          left_join(morphometry, by = "lake_id") %>%
                          left_join(landuse, by = "lake_id") %>%
                          left_join(soil, by = "lake_id") %>%
                          left_join(climate, by = "lake_id") %>%
                          left_join(waterchem, by = "lake_id"), "Environment (all)", arrows_taxa = 1.6))

(pca_env_plots <- wrap_plots(pca_geography, pca_morphometry, pca_landuse, pca_soil, pca_climate, pca_waterchem,
                             pca_env_all) + plot_layout(guides = "collect"))

# Define function to plot PCA with dark theme
plotPCADark <- function(env_data, plot_title, arrows_taxa = 1) {
  # Format environmental data
  env_data <- env_data %>%
    filter(lake_id %in% lakes) %>%
    dplyr::select(lake_id, where(is.numeric)) %>%
    column_to_rownames("lake_id")
  
  # Scale and center data
  env_sc <- scale(env_data, center = TRUE, scale = TRUE)
  
  # Compute PCA
  pca <- rda(env_sc)
  #summary(pca)  # Summarize PCA results
  
  # Extract PCA site and species scores
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "lake_id") %>%
    left_join(metadata, by = "lake_id")
  pca_species <- scores(pca)$species %>%
    as_tibble(rownames = "variable") %>%
    left_join(variable_labels_nounits, by = "variable")
  
  # Calculate variation explained by PCA dimensions 1 and 2
  pc1_pctvar <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1)
  pc2_pctvar <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1)
  
  # Plot PCA
  pca_plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#E6E6E6") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#E6E6E6") +
    geom_point(aes(x = PC1, y = PC2,
                   colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic")),
                   shape = factor(ecozone, levels = ecozone_longitude_order)),
               data = pca_sites, size = 5, alpha = 0.9, stroke = 2) +
    geom_text(aes(x = (arrows_taxa + 0.1)*PC1, y = (arrows_taxa + 0.1)*PC2, label = variable_curated),
              data = pca_species,
              colour = "white", size = 4) +
    geom_segment(aes(x = 0, y = 0, xend = arrows_taxa*PC1, yend = arrows_taxa*PC2),
                 data = pca_species,
                 arrow = arrow(length = unit(x = 0.2, units = "cm")),
                 alpha = 0.7, colour = "#A0A0A0") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    scale_colour_manual(values = palette_trophic) +
    labs(x = paste("PC axis 1 (", pc1_pctvar, "%)", sep = ""),
         y = paste("PC axis 2 (", pc2_pctvar, "%)", sep = ""),
         colour = "Trophic state") +
    #ggtitle(plot_title) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "#383838"),
          legend.key = element_blank(), legend.position = "right",
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  return(pca_plot)
}

(pca_env_all_dark <- plotPCADark(geography %>%
                          dplyr::select(-contains("_x"), -contains("_y"), -longitude) %>%
                          left_join(morphometry, by = "lake_id") %>%
                          left_join(landuse, by = "lake_id") %>%
                          left_join(soil, by = "lake_id") %>%
                          left_join(climate, by = "lake_id") %>%
                          left_join(waterchem, by = "lake_id"), "Environment (all)", arrows_taxa = 1.6))
#ggsave("figures/lp2017-2019watercolumn18s_env_pca.pdf", pca_env_all_dark, "pdf", width = 12, height = 11, units = "in")
