# Plot sampling sites on map of Canada with ecozones

# Load libraries
library(rgdal)
library(broom)
library(tidyverse)
library(ggnewscale)

# Load palettes
source("scripts/00_palettes.R")


#### Import data ####
# Import metadata
source("scripts/05_metadata.R")


#### Format map data ####
# Import and format ecozone map data
ecozones <- readOGR(dsn = "data/ecozones", layer = "Canada_Ecozones_V5b_31s_15M_simplify")
proj4string(ecozones)
ecozones@data$id <- rownames(ecozones@data)
ecozone_dat <- as.data.frame(ecozones)
ecozones$id <- rownames(ecozone_dat)
ecozones_pts <- tidy(ecozones)
# Retain only the ecozones in dataset
ecozones_df <- merge(ecozones_pts, ecozones, by = "id", type = "left") %>%
  as_tibble() %>%
  filter(Name %in% unique(metadata$ecozone))

# Format Canada base map data
canada <- map_data(map = "world", region = c("Canada"))
coordinates(canada) <- ~long + lat
proj4string(canada) <- CRS("+proj=longlat +datum=NAD83")
canada <- spTransform(canada, CRS(proj4string(ecozones)))
identical(proj4string(canada),proj4string(ecozones))  # Should evaluate to TRUE
canada_map <- data.frame(canada)
names(canada_map)[names(canada_map) == "longitude"] <- "x"
names(canada_map)[names(canada_map) == "latitude"] <- "y"

# Format lake metadata for mapping
metadata_map <- metadata
coordinates(metadata_map) <- ~longitude + latitude
class(metadata_map)
proj4string(metadata_map) <- CRS("+proj=longlat +datum=NAD83")
metadata_map <- spTransform(metadata_map, CRS(proj4string(ecozones)))
identical(proj4string(metadata_map),proj4string(ecozones))  # Should evaluate to TRUE
metadata_map <- data.frame(metadata_map)
names(metadata_map)[names(metadata_map) == "longitude"] <- "x"
names(metadata_map)[names(metadata_map) == "latitude"] <- "y"


#### Plot map ####
(map_plot <- ggplot() +
    geom_polygon(data = canada_map, mapping = aes(x = long, y = lat, group = group),
                 fill = "#f2f2f2", size = 0.1) +
    geom_polygon(data = ecozones_df, aes(x = long, y = lat, group = group,
                                         fill = factor(Name, levels = ecozone_longitude_order))) +
   scale_fill_manual(values = palette_ecozone_spectral, name = "Ecozone") +
   new_scale("fill") +
   geom_point(aes(x = x, y = y,
                   fill = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic"))),
               data = metadata_map,
               colour = "black", alpha = 0.9, size = 2.6, pch = 21, stroke = 0.01) +
   scale_fill_manual(values = palette_trophic, name = "Trophic state") +
   theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
          panel.border = element_blank(),
          legend.position = "right") +
    coord_equal(ratio = 1, xlim = c(-2200000, 2800000), ylim = c(-600000, 3700000)))  # Wide map limits
#ggsave("figures/lp2017-2019watercolumn18s_map.pdf", map_plot, "pdf", width = 10, height = 10, units = "in")
