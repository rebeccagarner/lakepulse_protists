# Plot diversity and environmental data radial bar plots

# Load libraries
library(tidyverse)
library(ggnewscale)
library(dendroextras)
library(ggpubr)

# Load palettes
source("scripts/00_palettes.R")


#### Import data ####
# Import metadata
source("scripts/05_metadata.R")

# Import melted protist ASV table
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_function.tsv", col_names = TRUE)
length(unique(asv_melt$lake_id))  # Number of lakes


#### Format environmental data ####
# Select land use fraction variables
landuse_fractions <- landuse %>%
  select(lake_id,
         fraction_agriculture, fraction_pasture, fraction_forestry,
         fraction_built, fraction_natlandscapes)

# Cluster lakes by land use fractions
landuse_fractions_df <- landuse_fractions %>%
  column_to_rownames("lake_id")

landuse_fractions_dist <- dist(landuse_fractions_df, method = "euclidean")
landuse_fractions_clust <- hclust(landuse_fractions_dist, method = "ward.D2")
landuse_fractions_dendro <- as.dendrogram(landuse_fractions_clust)
#plot(landuse_fractions_dendro)

landuse_fractions_order <- tibble(lake_id = labels(landuse_fractions_dendro)) %>%
  mutate(lake_order = row_number()) %>%
  mutate(lake_order = as.character(str_pad(lake_order, max(nchar(lake_order)), "left", 0)))

# Assign heatmap y-coordinates to land use variables
# Join cluster order
landuse_fractions <- landuse_fractions %>%
   pivot_longer(!lake_id, names_to = "landuse", values_to = "pct_landuse") %>%
   mutate(plot_index = case_when(landuse == "fraction_natlandscapes" ~ 5,
                                 landuse == "fraction_forestry" ~ 5.5,
                                 landuse == "fraction_built" ~ 6,
                                 landuse == "fraction_pasture" ~ 6.5,
                                 landuse == "fraction_agriculture" ~ 7)) %>%
   left_join(landuse_fractions_order, by = "lake_id")

# Join cluster order to metadata
metadata_order <- metadata %>%
   left_join(landuse_fractions_order, by = "lake_id")


#### Format ASV data ####
# Calculate relative abundance of divisions by lake
sum(asv_melt$relseqs)

divisions_relseqs <- asv_melt %>%
   mutate(division = case_when(!is.na(supergroup) & is.na(division) ~ str_c("Unclassified ", supergroup),
                               TRUE ~ division)) %>%
   group_by(lake_id, division) %>%
   summarize(pct_relseqs = sum(relseqs) * 100) %>%
   left_join(landuse_fractions_order, by = "lake_id")

# Calculate relative abundance of functional groups by lake
functions_relseqs <- asv_melt %>%
   group_by(lake_id, trophic_group) %>%
   summarize(pct_relseqs = sum(relseqs) * 100) %>%
   left_join(landuse_fractions_order, by = "lake_id")


#### Plot heatmap ####
# Plot land use heatmap
(landuse_heatmap <- ggplot() +
    geom_tile(aes(x = lake_order, y = plot_index, fill = pct_landuse),
              data = landuse_fractions %>%
                 filter(landuse == "fraction_natlandscapes"),
              colour = "#f6f6f6", height = 0.5) +
    scale_fill_gradient(low = "white", high = "#469E6F", name = "Natural landscapes (%)") +
    new_scale("fill") +
    geom_tile(aes(x = lake_order, y = plot_index, fill = pct_landuse),
              data = landuse_fractions %>%
                 filter(landuse == "fraction_forestry"),
              colour = "#f6f6f6", height = 0.5) +
    scale_fill_gradient(low = "white", high = "#CB6C43", name = "Forestry (%)") +
    new_scale("fill") +
    geom_tile(aes(x = lake_order, y = plot_index, fill = pct_landuse),
              data = landuse_fractions %>%
                 filter(landuse == "fraction_built"),
              colour = "#f6f6f6", height = 0.5) +
    scale_fill_gradient(low = "white", high = "#E63946", "Built (%)") +
    new_scale("fill") +
    geom_tile(aes(x = lake_order, y = plot_index, fill = pct_landuse),
              data = landuse_fractions %>%
                 filter(landuse == "fraction_pasture"),
              colour = "#f6f6f6", height = 0.5) +
    scale_fill_gradient(low = "white", high = "#FADD00", "Pasture (%)") +
    new_scale("fill") +
    geom_tile(aes(x = lake_order, y = plot_index, fill = pct_landuse),
              data = landuse_fractions %>%
                 filter(landuse == "fraction_agriculture"),
              colour = "#f6f6f6", height = 0.5) +
    scale_fill_gradient(low = "white", high = "#F28218", "Agriculture (%)") +
    coord_polar(theta = "x", start = 3.927, direction = -1) +
    scale_y_continuous(limits = c(-30, max(landuse_fractions$plot_index) + 1)) +  # Decrease x limit to narrow the circle width
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()))
#ggsave("figures/lp2017-2019watercolumn18s_radialheatmap_landuse.pdf", landuse_heatmap, device = "pdf", width = 18, height = 18, units = "in")


#### Plot bar plots ####
# Plot taxonomy bar plots ----
(taxonomy_radialbarplot <- ggplot() +      
   geom_bar(aes(x = lake_order, y = pct_relseqs,
                fill = factor(division, levels = c("Ochrophyta", "Opalozoa", "Pseudofungi", "Sagenista", "Stramenopiles_X",
                                                   #"Metazoa", "Fungi"
                                                   "Cryptophyta", "Katablepharidophyta", "Centroheliozoa", "Haptophyta", "Telonemia",
                                                   "Dinoflagellata", "Ciliophora", "Apicomplexa", "Perkinsea", "Protalveolata_X",
                                                   "Chlorophyta", "Streptophyta", "Rhodophyta", "Glaucophyta",
                                                   "Cercozoa", "Foraminifera",
                                                   "Choanoflagellida", "Mesomycetozoa",  # Opisthokonts
                                                   "Discoba", "Metamonada", "Malawimonadidae",
                                                   "Lobosa", "Conosa", "Breviatea",
                                                   "Hilomonadea", "Apusomonadidae"))),
            data = divisions_relseqs, stat = "identity") +
   scale_fill_manual(values = palette_division, na.value = "black") +
   coord_polar(theta = "x", start = 3.927, direction = -1) +
   scale_y_continuous(limits = c(-300,101)) +  # Decrease x limit to narrow the circle width
   theme(axis.text = element_blank(), axis.ticks = element_blank(),
         panel.background = element_blank(), panel.grid = element_blank(),
         axis.title = element_blank(),
         legend.position = "none", legend.key = element_blank()) +
   labs(fill = "Composition (% sequences)"))
#ggsave("figures/lp2017-2019watercolumn18s_radialbarplot_taxonomy.pdf", taxonomy_radialbarplot, device = "pdf", width = 16, height = 16, units = "in")

# Extract bar plot legend
# First need to change legend position to "right" in above plot
(radialbarplot_legend <- as_ggplot(get_legend(radialbarplot)))
#ggsave("figures/lp2017-2019watercolumn18s_radialbarplot_taxonomy_legend.pdf", radialbarplot_legend, device = "pdf", width = 3, height = 6, units = "in")

# Plot functional group bar plots ----
(function_radialbarplot <- ggplot() +      
   geom_bar(aes(x = lake_order, y = pct_relseqs,
                fill = factor(trophic_group, levels = c("phototrophs, consumers, mixotrophs", "phototrophs", "phototrophs, mixotrophs",
                                                        "mixotrophs", "consumers, mixotrophs", "consumers", "consumers, parasites",
                                                        "parasites"))),
            data = functions_relseqs, stat = "identity") +
   scale_fill_manual(values = palette_function, na.value = "#EEEEEE", name = "Functional group (% sequences)") +
   geom_point(aes(x = lake_order, y = 112,
                  colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
              data = metadata_order, size = 2.5) +
   scale_colour_manual(values = palette_trophic, name = "Trophic status") +
   coord_polar(theta = "x", start = 3.927, direction = -1) +
   scale_y_continuous(limits = c(-800,115)) +  # Decrease x limit to narrow the circle width
   theme(axis.text = element_blank(), axis.ticks = element_blank(),
         panel.background = element_blank(), panel.grid = element_blank(),
         axis.title = element_blank(),
         legend.position = "right", legend.key = element_blank()))
#ggsave("figures/lp2017-2019watercolumn18s_radialbarplot_function.pdf", function_radialbarplot, device = "pdf", width = 18, height = 18, units = "in")
