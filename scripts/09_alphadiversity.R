# Perform alpha-diversity analyses

# Load libraries
library(tidyverse)
library(vegan)
library(adespatial)
library(ape)
library(ggpubr)
library(ade4)
library(ggpmisc)
library(picante)
library(ggnewscale)
library(rgdal)
library(broom)
library(Hmisc)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import (generalized alpha = 0.5) UniFrac distance matrix
load("output/dissim/lp2017-2019watercolumn18s_protists_unifrac_gen0pt5alpha.rda")

# Import melted ASV data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)

# Import rooted phylogenetic tree
tree <- read.tree("output/trees/lp2017-2019watercolumn18s_protists_366lakes_sina_gtr_rooted.tree")
length(unique(asv_melt$asv_code)) == length(tree$tip.label)  # Should evaluate to TRUE

# Import metadata
source("scripts/05_metadata.R")


#### Create taxonomy table ####
taxonomy <- asv_melt %>%
  distinct(asv_code, kingdom, supergroup, division, class, order, family, genus, species)


#### Create ASV tables ####
# ASV table of relative sequence abundance
asvtable <- asv_melt %>%
  select(lake_id, asv_code, relseqs) %>%
  spread(asv_code, relseqs, fill = 0) %>%
  column_to_rownames("lake_id")

# ASV table of sequence counts
asvtable_nseqs <- asv_melt %>%
  select(lake_id, asv_code, nseqs) %>%
  spread(asv_code, nseqs, fill = 0) %>%
  column_to_rownames("lake_id")


#### Subsample ASV data to lowest sample sequence abundance ####
# Calculate minimum sequence abundance
(nseqs_min <- asv_melt %>%
   select(lake_id, asv_code, nseqs) %>%
   distinct(lake_id, asv_code, .keep_all = TRUE) %>%
   group_by(lake_id) %>%
   dplyr::summarize(total_nseqs = sum(nseqs)) %>%
   slice_min(total_nseqs) %>%
   pull(total_nseqs))

# Subsample ASVs to minimum sequence abundance  (only run once!)
#asvtable_rarefied <- rrarefy(asvtable_nseqs, nseqs_min)
#sum(asvtable_rarefied)/nseqs_min  # Should evaluate to the number of lakes

# Save rarefied ASV table
# asvtable_rarefied %>%
#   as_tibble(rownames = "lake_id") %>%
#   pivot_longer(!lake_id, names_to = "asv_code", values_to = "nseqs") %>%
#   filter(nseqs > 0) %>%
#   write_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_rarefied.tsv", col_names = TRUE)


#### Import and format rarefied ASV table ####
# Import rarefied ASV data
asvtable_rarefied <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_rarefied.tsv", col_names = TRUE) %>%
  pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0) %>%
  column_to_rownames("lake_id")
dim(asvtable_rarefied)


#### Calculate alpha-diversity indices ####
# Richness ----
richness <- specnumber(asvtable_rarefied)
range(richness)

# Shannon diversity index ----
shannon <- diversity(asvtable_rarefied, "shannon")
range(shannon)

# Pielou's evenness index ----
pielou <- shannon/log(richness)
range(pielou)

# Faith's phylogenetic diversity index ----
faith <- pd(asvtable_rarefied, tree, include.root = FALSE) %>%
  as_tibble(rownames = "lake_id") %>%
  select(-SR) %>%
  rename(faith = PD)
range(faith$faith)

# Summarize alpha-diversity indices
diversity <- tibble(lake_id = rownames(asvtable_rarefied), richness, shannon, pielou) %>%
  left_join(faith, by = "lake_id")

# Save alpha-diversity table to file
# diversity %>%
#   write_tsv("output/diversity/lp2017-2019watercolumn18s_alphadiversity.tsv", col_names = TRUE)

# Load saved alpha-diversity table
diversity <- read_tsv("output/diversity/lp2017-2019watercolumn18s_alphadiversity.tsv", col_names = TRUE)


#### Examine alpha-diversity ####
# Correlate alpha-diversity indices
plotCorr(diversity)

(richness_faith_plot <- diversity %>%
    left_join(metadata, by = "lake_id") %>%
    ggplot(aes(x = richness, y = faith)) +
    geom_point(aes(colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic")),
                   shape = factor(ecozone, levels = ecozone_longitude_order)),
               size = 3, alpha = 0.9, stroke = 2) +
    labs(x = "ASV richness",
         y = "Faith's phylogenetic diversity index") +
    scale_colour_manual(values = palette_trophic, name = "Trophic state") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    theme_bw())
#ggsave("figures/lp2017-2019watercolumn18s_richness_faithpd_plot.pdf", richness_faith_plot, "pdf", width = 10, height = 8, units = "in")

# Diversity along east-west gradient
(shannon_eastwest_plot <- diversity %>%
    left_join(metadata, by = "lake_id") %>%
    ggplot() +
    geom_point(aes(x = longitude, y = shannon,
                   shape = factor(ecozone, levels = ecozone_longitude_order),
                   colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
               size = 3, alpha = 0.9, stroke = 2) +
    labs(x = "Longitude (°W)",
         y = "Shannon's index") +
    scale_colour_manual(values = palette_trophic, name = "Trophic state") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    theme_bw())
#ggsave("figures/lp2017-2019watercolumn18s_shannon_eastwest_plot.pdf", shannon_eastwest_plot, "pdf", width = 10, height = 8, units = "in")

(faithpd_eastwest_plot <- diversity %>%
    left_join(metadata, by = "lake_id") %>%
    ggplot() +
    geom_point(aes(x = longitude, y = faith,
                   shape = factor(ecozone, levels = ecozone_longitude_order),
                   colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
               size = 3, alpha = 0.9, stroke = 2) +
    labs(x = "Longitude (°W)",
         y = "Faith's phylogenetic diversity index") +
    scale_colour_manual(values = palette_trophic, name = "Trophic state") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    theme_bw())
#ggsave("figures/lp2017-2019watercolumn18s_faith_eastwest_plot.pdf", faithpd_eastwest_plot, "pdf", width = 10, height = 8, units = "in")

# Diversity by trophic state
trophic <- metadata %>%
  dplyr::select(lake_id, trophic_status)

(diversity_trophic_boxplot <- diversity %>%
    pivot_longer(!lake_id, names_to = "diversity_index", values_to = "value") %>%
    left_join(trophic, by = "lake_id") %>%
    left_join(ecozones, by = "lake_id") %>%
    ggplot() +
    facet_wrap(~factor(diversity_index, levels = c("richness", "shannon", "pielou", "faith")),
               ncol = 2, scales = "free") +
    geom_boxplot(aes(x = factor(trophic_status, levels = c("ultraoligotrophic", "oligotrophic", "mesotrophic", "mesoeutrophic", "eutrophic", "hypereutrophic")),
                     y = value,
                     colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
                 outlier.shape = NA) +
    geom_jitter(aes(x = factor(trophic_status, levels = c("ultraoligotrophic", "oligotrophic", "mesotrophic", "mesoeutrophic", "eutrophic", "hypereutrophic")),
                    y = value,
                    colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic")),
                    shape = factor(ecozone, levels = ecozone_longitude_order)),
                height = 0, alpha = 0.7, size = 2) +
    coord_flip() +
    scale_colour_manual(values = palette_trophic, name = "Trophic state") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    theme_light() %+replace%
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.background = element_rect(fill = "black"),
          strip.text = element_text(face = "bold", colour = "white"),
          legend.position = "right"))
#ggsave("figures/lp2017-2019watercolumn18s_alphadiversity_trophic_boxplots.pdf", diversity_trophic_boxplot, "pdf", width = 14, height = 6, units = "in")

# Calculate mean alpha-diversity by trophic state
diversity %>%
  pivot_longer(!lake_id, names_to = "diversity_index", values_to = "value") %>%
  left_join(trophic, by = "lake_id") %>%
  group_by(trophic_status, diversity_index) %>%
  dplyr::summarize(mean_value = mean(value)) %>%
  arrange(-mean_value)

# Perform ANOVA to test whether richness is different across trophic states
richness_anova <- aov(richness ~ trophic_status, data = diversity %>%
                        left_join(trophic, by = "lake_id"))
summary(richness_anova)

tukey_richness <- TukeyHSD(richness_anova)$trophic_status %>%
  as_tibble(rownames = "comparison") %>%
  mutate(distance = "richness")

# Perform ANOVA to test whether Shannon diversity is different across trophic states
shannon_anova <- aov(shannon ~ trophic_status, data = diversity %>%
                       left_join(trophic, by = "lake_id"))
summary(shannon_anova)

tukey_shannon <- TukeyHSD(shannon_anova)$trophic_status %>%
  as_tibble(rownames = "comparison") %>%
  mutate(distance = "shannon")

# Perform ANOVA to test whether Pielou's evenness index is different across trophic states
pielou_anova <- aov(pielou ~ trophic_status, data = diversity %>%
                      left_join(trophic, by = "lake_id"))
summary(pielou_anova)

tukey_pielou <- TukeyHSD(pielou_anova)$trophic_status %>%
  as_tibble(rownames = "comparison") %>%
  mutate(distance = "pielou")

# Perform ANOVA to test whether Faith's diversity index is different across trophic states
faith_anova <- aov(faith ~ trophic_status, data = diversity %>%
                     left_join(trophic, by = "lake_id"))
summary(faith_anova)

tukey_faith <- TukeyHSD(faith_anova)$trophic_status %>%
  as_tibble(rownames = "comparison") %>%
  mutate(distance = "faith")

# Combine Tukey's test results
tukey_diversity_all <- bind_rows(tukey_richness,
                                 tukey_shannon,
                                 tukey_pielou,
                                 tukey_faith) %>%
  clean_names() %>%
  filter(!grepl("ultraoligotrophic", comparison)) %>%
  filter(p_adj < 0.05) %>%
  arrange(factor(distance, levels = c("richness", "shannon", "pielou", "faith")), p_adj) %>%
  mutate(significance_level = case_when(p_adj < 0.05 & p_adj >= 0.01 ~ "*",
                                        p_adj < 0.01 & p_adj >= 0.001 ~ "**",
                                        p_adj < 0.001 ~ "***")) %>%
  select(distance, comparison, diff, lwr, upr, p_adj, significance_level) %>%
  rename(difference_in_means = diff,
         lower_endpoint_interval = lwr,
         upper_endpoint_interval = upr)

# tukey_diversity_all %>%
#   write_csv("output/diversity/lp2017-2019watercolumn18s_alphadiversity_trophic_tukey.csv")

# Diversity by ecozones
ecozones <- metadata %>%
  dplyr::select(lake_id, ecozone)

(diversity_ecozone_boxplot <- diversity %>%
    pivot_longer(!lake_id, names_to = "diversity_index", values_to = "value") %>%
    left_join(ecozones, by = "lake_id") %>%
    ggplot() +
    facet_wrap(~factor(diversity_index, levels = c("richness", "shannon", "pielou", "faith")),
               ncol = 1, scales = "free") +
    geom_boxplot(aes(x = factor(ecozone, levels = ecozone_longitude_order),
                     y = value,
                     colour = factor(ecozone, levels = ecozone_longitude_order))) +
    geom_jitter(aes(x = factor(ecozone, levels = ecozone_longitude_order),
                    y = value,
                    colour = factor(ecozone, levels = ecozone_longitude_order)),
                alpha = 0.1) +
    scale_colour_manual(values = palette_ecozone_spectral, name = "Ecozone") +
    theme_light() %+replace%
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "black")))
#ggsave("figures/lp2017-2019watercolumn18s_alphadiversity_ecozone_boxplots.pdf", diversity_ecozone_boxplot, "pdf", width = 8, height = 15, units = "in")


#### Plot correlations between diversity metrics and environmental variables #####
# Define function to plot correlations with select data
plotCorrSelect <- function(data, metadata) {
  data_metadata <- data %>%
    left_join(metadata, by = "lake_id") %>%
    select(-latitude, -longitude)
  
  cor(data_metadata[,-1]) %>%
    as_tibble(rownames = "var1") %>%
    pivot_longer(!var1, names_to = "var2", values_to = "correlation") %>%
    filter(var1 %in% names(data)) %>%
    filter(!var2 %in% names(diversity)) %>%
    ggplot() +
    geom_tile(aes(x = var1, y = reorder(var2, correlation), fill = correlation)) +
    geom_text(aes(x = var1, y = reorder(var2, correlation), label = round(correlation, 2))) +
    scale_fill_gradient2(low = "red", mid = "white", high = "dodgerblue", midpoint = 0, limits = c(-1,1)) +
    scale_x_discrete(position = "top") +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0),
          panel.background = element_blank(), axis.ticks = element_blank())
}

# Define function to add lat-long coordinates to environmental data matrix
addLatLong <- function(envdata) {
  env_latlong <- latlong %>%
    left_join(envdata, by = "lake_id")
  
  return(env_latlong)
}

# Add lat-long coordinate information to environmental matrices (except geography)
morphometry_latlong <- addLatLong(morphometry)
landuse_latlong <- addLatLong(landuse)
soil_latlong <- addLatLong(soil)
landuse_soil <- landuse %>%
  left_join(soil, by = "lake_id")
landuse_soil_latlong <- addLatLong(landuse_soil)
climate_latlong <- addLatLong(climate)
waterchem_latlong <- addLatLong(waterchem)

# Plot correlations and compute correlation p-values
plotCorrSelect(diversity, geography)
rcorr(as.matrix(diversity[,-1]), as.matrix(geography[,-1]))$P
plotCorrSelect(diversity, morphometry_latlong)
rcorr(as.matrix(diversity[,-1]), as.matrix(morphometry_latlong[,-1]))$P
plotCorrSelect(diversity, landuse_latlong)
rcorr(as.matrix(diversity[,-1]), as.matrix(landuse_latlong[,-1]))$P
plotCorrSelect(diversity, soil_latlong)
rcorr(as.matrix(diversity[,-1]), as.matrix(soil_latlong[,-1]))$P
plotCorrSelect(diversity, landuse_soil_latlong)
rcorr(as.matrix(diversity[,-1]), as.matrix(landuse_soil_latlong[,-1]))$P
plotCorrSelect(diversity, climate_latlong)
rcorr(as.matrix(diversity[,-1]), as.matrix(climate_latlong[,-1]))$P
plotCorrSelect(diversity, waterchem_latlong)
rcorr(as.matrix(diversity[,-1]), as.matrix(waterchem_latlong[,-1]))$P

# Plot diversity metrics vs. correlated environmental variables
diversity_env <- diversity %>%
  left_join(env_all, by = "lake_id") %>%
  select(-latitude, -longitude, -contains("mem"))

vars_cor <- cor(diversity_env[,-1]) %>%
  as_tibble(rownames = "var1") %>%
  pivot_longer(!var1, names_to = "var2", values_to = "correlation") %>%
  filter(var1 %in% names(diversity)) %>%
  filter(!var2 %in% names(diversity)) %>%
  filter(abs(correlation) >= 0.1) %>%
  distinct(var2) %>%
  pull(var2)

(diversity_env_lmplots <- diversity_env %>%
    pivot_longer(!names(diversity), names_to = "var2", values_to = "value2") %>%
    pivot_longer(!c(lake_id, var2, value2), names_to = "var1", values_to = "value1") %>%
    left_join(ecozones, by = "lake_id") %>%
    left_join(trophic, by = "lake_id") %>%
    ggplot(aes(x = value2, y = value1)) +
    facet_grid(var1~var2, scales = "free") +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, colour = "lightblue") +
    stat_poly_eq(formula = y ~ x,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE, label.y = "bottom") +
    geom_point(aes(colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic")),
                   shape = factor(ecozone, levels = ecozone_longitude_order)),
               size = 1, alpha = 0.5, stroke = 2) +
    scale_colour_manual(values = palette_trophic, name = "Trophic state") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    theme_bw() %+replace%
    theme(axis.title = element_blank()))
#ggsave("figures/lp2017-2019watercolumn18s_diversity_env_lmscatterplots.pdf", diversity_env_lmplots, "pdf", width = 140, height = 20, units = "in", limitsize = FALSE)


#### Map alpha-diversity ####
# Format Canada/USA base map data
ecozones <- readOGR(dsn = "data/ecozones", layer = "Canada_Ecozones_V5b_31s_15M_simplify")
proj4string(ecozones)
ecozones@data$id <- rownames(ecozones@data)
ecozone_dat <- as.data.frame(ecozones)
ecozones$id <- rownames(ecozone_dat)
ecozones_pts <- tidy(ecozones)
ecozones_df <- merge(ecozones_pts, ecozones, by = "id", type = "left") %>%
  as_tibble() %>%
  filter(Name %in% unique(metadata$ecozone))

canada <- map_data(map = "world", region = c("Canada"))
coordinates(canada) <- ~long + lat
proj4string(canada) <- CRS("+proj=longlat +datum=NAD83")
canada <- spTransform(canada, CRS(proj4string(ecozones)))
identical(proj4string(canada),proj4string(ecozones))  # Should evaluate to TRUE
canada_map <- data.frame(canada)
names(canada_map)[names(canada_map) == "longitude"] <- "x"
names(canada_map)[names(canada_map) == "latitude"] <- "y"

coordinates(diversity_env) <- ~longitude_x + latitude_y
class(diversity_env)
proj4string(diversity_env) <- CRS("+proj=longlat +datum=NAD83")
diversity_env <- spTransform(diversity_env, CRS(proj4string(ecozones)))
identical(proj4string(diversity_env),proj4string(ecozones))  # Should evaluate to TRUE
diversity_env <- data.frame(diversity_env)
names(diversity_env)[names(diversity_env) == "longitude_x"] <- "x"
names(diversity_env)[names(diversity_env) == "latitude_y"] <- "y"

richness_jenks <- classInt::classIntervals(diversity_env$richness, n = 5, style = "jenks")$brks %>%
  round(1)

diversity_env <- diversity_env %>%
  mutate(richness_order = case_when(richness <= richness_jenks[2] ~ 1,
                                    richness > richness_jenks[2] & richness <= richness_jenks[3] ~ 2,
                                    richness > richness_jenks[3] & richness <= richness_jenks[4] ~ 3,
                                    richness > richness_jenks[4] & richness <= richness_jenks[5] ~ 4,
                                    richness > richness_jenks[5] & richness <= richness_jenks[6] ~ 5)) %>%
  mutate(richness_jenks = case_when(richness <= richness_jenks[2] ~ paste0(richness_jenks[1], " - ", richness_jenks[2]),
                                    richness > richness_jenks[2] & richness <= richness_jenks[3] ~ paste0(richness_jenks[2], " - ", richness_jenks[3]),
                                    richness > richness_jenks[3] & richness <= richness_jenks[4] ~ paste0(richness_jenks[3], " - ", richness_jenks[4]),
                                    richness > richness_jenks[4] & richness <= richness_jenks[5] ~ paste0(richness_jenks[4], " - ", richness_jenks[5]),
                                    richness > richness_jenks[5] & richness <= richness_jenks[6] ~ paste0(richness_jenks[5], " - ", richness_jenks[6])))

(richness_map_plot <- ggplot() +
    geom_polygon(data = canada_map, mapping = aes(x = long, y = lat, group = group),
                 fill = "#f2f2f2", size = 0.1) +
    geom_polygon(data = ecozones_df, aes(x = long, y = lat, group = group, fill = Name),
                 alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = palette_ecozone_spectral) +
    new_scale("fill") +
    geom_point(aes(x = x, y = y, fill = fct_reorder(richness_jenks, richness_order)), data = diversity_env,
               alpha = 0.9, size = 2.6, pch = 21, stroke = 0.01) +
    scale_fill_manual(values = c("#4DAC26", "#B8E186", "#F7F7F7", "#F1B6DA", "#D01C8B"),
                      name = "Richness", guide = guide_legend(reverse = TRUE)) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
          panel.border = element_blank(),
          legend.position = "right", legend.key = element_blank()) +
    coord_equal(ratio = 1, xlim = c(-2200000, 2800000), ylim = c(-600000, 3700000)))
#ggsave("figures/lp2017-2019watercolumn18s_richness_map.pdf", richness_map_plot, "pdf", width = 10, height = 10, units = "in")

shannon_jenks <- classInt::classIntervals(diversity_env$shannon, n = 5, style = "jenks")$brks %>%
  round(1)

diversity_env <- diversity_env %>%
  mutate(shannon_order = case_when(shannon <= shannon_jenks[2] ~ 1,
                                   shannon > shannon_jenks[2] & shannon <= shannon_jenks[3] ~ 2,
                                   shannon > shannon_jenks[3] & shannon <= shannon_jenks[4] ~ 3,
                                   shannon > shannon_jenks[4] & shannon <= shannon_jenks[5] ~ 4,
                                   shannon > shannon_jenks[5] & shannon <= shannon_jenks[6] ~ 5)) %>%
  mutate(shannon_jenks = case_when(shannon <= shannon_jenks[2] ~ paste0(shannon_jenks[1], " - ", shannon_jenks[2]),
                                   shannon > shannon_jenks[2] & shannon <= shannon_jenks[3] ~ paste0(shannon_jenks[2], " - ", shannon_jenks[3]),
                                   shannon > shannon_jenks[3] & shannon <= shannon_jenks[4] ~ paste0(shannon_jenks[3], " - ", shannon_jenks[4]),
                                   shannon > shannon_jenks[4] & shannon <= shannon_jenks[5] ~ paste0(shannon_jenks[4], " - ", shannon_jenks[5]),
                                   shannon > shannon_jenks[5] & shannon <= shannon_jenks[6] ~ paste0(shannon_jenks[5], " - ", shannon_jenks[6])))

(shannon_map_plot <- ggplot() +
    geom_polygon(data = canada_map, mapping = aes(x = long, y = lat, group = group),
                 fill = "#f2f2f2", size = 0.1) +
    geom_polygon(data = ecozones_df, aes(x = long, y = lat, group = group, fill = Name),
                 alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = palette_ecozone_spectral) +
    new_scale("fill") +
    geom_point(aes(x = x, y = y, fill = fct_reorder(shannon_jenks, shannon_order)), data = diversity_env,
               alpha = 0.9, size = 2.6, pch = 21, stroke = 0.01) +
    scale_fill_manual(values = c("#4DAC26", "#B8E186", "#F7F7F7", "#F1B6DA", "#D01C8B"),
                      name = "Shannon index", guide = guide_legend(reverse = TRUE)) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
          panel.border = element_blank(),
          legend.position = "right", legend.key = element_blank()) +
    #scale_size_continuous(range = c(1,10)) +
    coord_equal(ratio = 1, xlim = c(-2200000, 2800000), ylim = c(-600000, 3700000)))
#ggsave("figures/lp2017-2019watercolumn18s_shannon_map.pdf", shannon_map_plot, "pdf", width = 10, height = 10, units = "in")
