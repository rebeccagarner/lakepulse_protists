# Plot environmental data distributions

# Load libraries
library(tidyverse)
library(patchwork)


#### Import and format data ####
# Import metadata
source("scripts/05_metadata.R")


#### Plot histograms ####
(geography_histogram <- geography %>%
    dplyr::select(lake_id, latitude, longitude, altitude) %>%
    plotDistrib(plot_ncols = 16))

(morphometry_histogram <- morphometry %>%
    plotDistrib(plot_ncols = 16))

(watershed_histogram <- landuse %>%
    left_join(soil, by = "lake_id") %>%
    plotDistrib(plot_ncols = 16))

(climate_histogram <- climate %>%
    plotDistrib(plot_ncols = 16))

(waterchem_histogram <- waterchem %>%
    plotDistrib(plot_ncols = 16))

# Combine histograms
histograms_env_all_layout <- "
AAA#############
BBBBBBBBB#######
CCCCCCCCCCCCCCCC
DDDDD###########
EEEEEEEEEEEEEE##
"
(histograms_env_plots_all <- geography_histogram /
    morphometry_histogram /
    watershed_histogram /
    climate_histogram /
    waterchem_histogram + 
    plot_layout(design = histograms_env_all_layout, guides = "collect"))
#ggsave("figures/lp2017-2019watercolumn18s_env_histograms_all.pdf", histograms_env_plots_all, "pdf", width = 32, height = 10, units = "in")
