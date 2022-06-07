# Evaluate local contributions to beta-diversity (LCBD)

# Load libraries
library(tidyverse)
library(vegan)
library(adespatial)
library(ape)
library(ggpubr)
library(ade4)
library(ggpmisc)
library(SoDA)
library(sp)
library(rgdal)
library(broom)
library(ggnewscale)
library(Hmisc)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import melted protists data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)
length(unique(asv_melt$asv_code))

# Import rarefied ASV data
asv_melt_rarefied <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_rarefied.tsv", col_names = TRUE)

# Import alpha-diversity metrics
diversity <- read_tsv("output/diversity/lp2017-2019watercolumn18s_alphadiversity.tsv", col_names = TRUE)

# Import metadata
source("scripts/05_metadata.R")


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

# Rarefied ASV table
asvtable_rarefied <- asv_melt_rarefied %>%
  pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0) %>%
  column_to_rownames("lake_id")


#### Local Contributions to Beta-Diversity ####
# Calculate beta-diversity based on Hellinger-transformed ASV compositions
betadiv <- beta.div(asvtable, method = "hellinger", nperm = 100)

# Assess total beta-diversity
betadiv$beta  # SStotal and BDtotal
lcbd <- tibble(lake_id = names(betadiv$LCBD),
               lcbd = betadiv$LCBD)

# Identify sites with significant Holm-corrected LCBD value
row.names(asvtable[which(p.adjust(betadiv$p.LCBD, "holm") <= 0.05),])

# Write LCBD to file
# lcbd %>%
#   write_tsv("output/diversity/lp2017-2019watercolumn18s_lcbd.tsv")

# Import LCBD
lcbd <- read_tsv("output/diversity/lp2017-2019watercolumn18s_lcbd.tsv", col_names = TRUE)


#### Correlate LCBD and environmental variables ####
# Define function to plot correlations with select data
plotCorrSelect <- function(data, metadata) {
  data_metadata <- data %>%
    left_join(metadata, by = "lake_id") %>%
    select(-latitude, -longitude)
  
  cor(data_metadata[,-1]) %>%
    as_tibble(rownames = "var1") %>%
    pivot_longer(!var1, names_to = "var2", values_to = "correlation") %>%
    filter(var1 %in% names(data)) %>%
    filter(!var2 %in% names(data)) %>%
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
plotCorrSelect(lcbd, geography)
rcorr(as.matrix(lcbd[,-1]), as.matrix(geography[,-1]))$P
plotCorrSelect(lcbd, morphometry_latlong)
rcorr(as.matrix(lcbd[,-1]), as.matrix(morphometry_latlong[,-1]))$P
plotCorrSelect(lcbd, landuse_latlong)
rcorr(as.matrix(lcbd[,-1]), as.matrix(landuse_latlong[,-1]))$P
plotCorrSelect(lcbd, soil_latlong)
rcorr(as.matrix(lcbd[,-1]), as.matrix(soil_latlong[,-1]))$P
plotCorrSelect(lcbd, landuse_soil_latlong)
rcorr(as.matrix(lcbd[,-1]), as.matrix(landuse_soil_latlong[,-1]))$P
plotCorrSelect(lcbd, climate_latlong)
rcorr(as.matrix(lcbd[,-1]), as.matrix(climate_latlong[,-1]))$P
plotCorrSelect(lcbd, waterchem_latlong)
rcorr(as.matrix(lcbd[,-1]), as.matrix(waterchem_latlong[,-1]))$P

# Plot diversity metrics vs. correlated environmental variables
lcbd_env <- lcbd %>%
  left_join(env_all, by = "lake_id") %>%
  select(-latitude, -longitude, -contains("mem"))

vars_cor <- cor(lcbd_env[,-1]) %>%
  as_tibble(rownames = "var1") %>%
  pivot_longer(!var1, names_to = "var2", values_to = "correlation") %>%
  filter(var1 %in% names(lcbd)) %>%
  filter(!var2 %in% names(lcbd)) %>%
  filter(abs(correlation) >= 0.1) %>%
  distinct(var2) %>%
  pull(var2)

(lcbd_env_lmplots <- lcbd_env %>%
    pivot_longer(!names(lcbd), names_to = "var2", values_to = "value2") %>%
    pivot_longer(!c(lake_id, var2, value2), names_to = "var1", values_to = "value1") %>%
    ggplot(aes(x = value2, y = value1)) +
    facet_grid(var1~var2, scales = "free") +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, colour = "lightblue") +
    stat_poly_eq(formula = y ~ x,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE, label.y = "bottom") +
    geom_point(alpha = 0.5) +
    theme_bw() %+replace%
    theme(axis.title = element_blank()))
#ggsave("figures/lp2017-2019watercolumn18s_lcbd_env_lmscatterplots.pdf", lcbd_env_lmplots, "pdf", width = 140, height = 5, units = "in", limitsize = FALSE)


#### Map LCBD ####
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

coordinates(lcbd_env) <- ~longitude_x + latitude_y
class(lcbd_env)
proj4string(lcbd_env) <- CRS("+proj=longlat +datum=NAD83")
lcbd_env <- spTransform(lcbd_env, CRS(proj4string(ecozones)))
identical(proj4string(lcbd_env),proj4string(ecozones))  # Should evaluate to TRUE
lcbd_env <- data.frame(lcbd_env)
names(lcbd_env)[names(lcbd_env) == "longitude_x"] <- "x"
names(lcbd_env)[names(lcbd_env) == "latitude_y"] <- "y"

lcbd_jenks <- classInt::classIntervals(lcbd_env$lcbd, n = 4, style = "jenks")$brks %>%
  round(4)

lcbd_env <- lcbd_env %>%
  mutate(lcbd_jenks = case_when(lcbd <= lcbd_jenks[2] ~ paste0(lcbd_jenks[1], " - ", lcbd_jenks[2]),
                                lcbd > lcbd_jenks[2] & lcbd <= lcbd_jenks[3] ~ paste0(lcbd_jenks[2], " - ", lcbd_jenks[3]),
                                lcbd > lcbd_jenks[3] & lcbd <= lcbd_jenks[4] ~ paste0(lcbd_jenks[3], " - ", lcbd_jenks[4]),
                                lcbd > lcbd_jenks[4] & lcbd <= lcbd_jenks[5] ~ paste0(lcbd_jenks[4], " - ", lcbd_jenks[5])))

(lcbd_map_plot <- ggplot() +
    geom_polygon(data = canada_map, mapping = aes(x = long, y = lat, group = group),
                 fill = "#f2f2f2", size = 0.1) +
    geom_polygon(data = ecozones_df, aes(x = long, y = lat, group = group, fill = Name),
                 alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = palette_ecozone_spectral) +
    new_scale("fill") +
    geom_point(aes(x = x, y = y, fill = lcbd_jenks), data = lcbd_env,
               alpha = 0.9, size = 2.6, pch = 21, stroke = 0.01) +
    scale_fill_manual(values = c("#4DAC26", "#B8E186", "#F1B6DA", "#D01C8B"),
                      name = "LCBD", guide = guide_legend(reverse = TRUE)) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
          panel.border = element_blank(),
          legend.position = "right", legend.key = element_blank()) +
    #scale_size_continuous(range = c(1,10)) +
    coord_equal(ratio = 1, xlim = c(-2200000, 2800000), ylim = c(-600000, 3700000)))
#ggsave("figures/lp2017-2019watercolumn18s_lcbd_hellinger_map.pdf", lcbd_map_plot, "pdf", width = 10, height = 10, units = "in")


#### Alpha-diversity indices ####
# Join alpha-diversity and LCBD
diversity <- diversity %>%
  left_join(lcbd, by = "lake_id")

# diversity %>%
#   left_join(metadata, by = "lake_id") %>%
#   write_tsv("output/diversity/lp2017-2019watercolumn18s_alphadiversity_lcbd.tsv", col_names = TRUE)

# Correlate LCBD and alpha-diversity indices
plotCorr(diversity)

(lcbd_richness_plot <- diversity %>%
    left_join(metadata, by = "lake_id") %>%
    ggplot(aes(x = richness, y = lcbd)) +
    geom_point(aes(shape = factor(ecozone, levels = ecozone_longitude_order),
                   colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
               size = 3, alpha = 0.9, stroke = 2) +
    labs(x = "ASV richness",
         y = "LCBD") +
    scale_colour_manual(values = palette_trophic, name = "Trophic status") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    theme_bw())
#ggsave("figures/lp2017-2019watercolumn18s_lcbd_richness_plot.pdf", lcbd_richness_plot, "pdf", width = 10, height = 8, units = "in")


#### Decompose beta-diversity ####
# Compute beta-diversity as percentage differences (Sorensen index + quantitative)
betadiv_bray <- beta.div.comp(asvtable_nseqs, coef = "S", quant = TRUE)

# Decompose beta-diversity into replacement (Repl) and abundance difference (RichDif)
betadiv_bray$part  # Beta-diversity partitioning vector

betadiv_total <- betadiv_bray$D  # Total beta-diversity
betadiv_replacement <- betadiv_bray$repl  # Replacement component of beta-diversity
betadiv_richnessdiff <- betadiv_bray$rich  # Richness difference component of beta-diversity

# Define function to plot PCoA
plotPCoA <- function(dist, plot_title) {
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
               data = pcoa_sites, size = 4, alpha = 0.9, stroke = 2) +
    scale_colour_manual(values = palette_trophic, name = "Trophic status") +
    scale_shape_manual(values = palette_ecozone_letters, name = "Ecozone") +
    labs(x = paste("PCoA axis 1 (", pcoa1_pctvar, "%)", sep = ""),
         y = paste("PCoA axis 2 (", pcoa2_pctvar, "%)", sep = "")) +
    ggtitle(plot_title) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          legend.key = element_blank(), legend.position = "right",
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  return(pcoa_plot)
}

# Transform beta-diversity components into rectangular matrices by PCoA
# PCoA of total beta-diversity ----
is.euclid(dist(betadiv_total))
(pcoa_betadiv_total <- plotPCoA(betadiv_total, "Total beta-diversity"))

# PCoA of replacement component of beta-diversity ----
is.euclid(dist(betadiv_replacement))
(pcoa_betadiv_replacement <- plotPCoA(betadiv_replacement, "Replacement"))

# PCoA of richness difference of beta-diversity ----
is.euclid(dist(betadiv_richnessdiff))
(pcoa_betadiv_richnessdiff <- plotPCoA(betadiv_richnessdiff, "Richness difference"))

# Combine PCoA plots of beta-diversity components
(pcoa_betadiv_plots <- ggarrange(pcoa_betadiv_total, pcoa_betadiv_replacement, pcoa_betadiv_richnessdiff,
                                 ncol = 3, nrow = 1,
                                 common.legend = TRUE, legend = "right"))
#ggsave("figures/lp2017-2019watercolumn18s_protists_betadiversity_partitioning_pcoa.pdf", pcoa_betadiv_plots, "pdf", width = 25, height = 8, units = "in")
