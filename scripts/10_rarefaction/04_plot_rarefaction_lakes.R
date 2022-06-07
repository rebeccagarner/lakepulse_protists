# Plot ASV rarefaction curves (all lakes)

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import melted protists data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)

# Import rarefaction tables
rarefaction_path <- "output/diversity/rarefaction/"
(rarefaction_files <- dir(path = rarefaction_path, pattern = "*_rarefaction.tsv"))
rarefaction <- tibble(file_name = rarefaction_files) %>%
  mutate(lake_id = str_remove(file_name, "lp2017-2019watercolumn18s_protists_")) %>%
  mutate(lake_id = str_remove(lake_id, "_rarefaction.tsv")) %>%
  mutate(file_contents = map(file_name, ~ read_tsv(file.path(rarefaction_path, .), col_names = TRUE))) %>%
  unnest(c(file_contents)) %>%
  select(-file_name)

# Import metadata
source("scripts/05_metadata.R")


#### Plot rarefaction curves ####
# Join lake ecozone and trophic status information to rarefaction table
lakes_ecozones_trophic <- metadata %>%
  dplyr::select(lake_id, ecozone, trophic_status)

rarefaction <- rarefaction %>%
  left_join(lakes_ecozones_trophic, by = "lake_id")

# Plot rarefaction curves
(rarefaction_curve <- rarefaction %>%
    ggplot(aes(x = rarefied, y = richness,
               colour = factor(trophic_status, levels = c("hypereutrophic", "eutrophic", "mesoeutrophic", "mesotrophic", "oligotrophic", "ultraoligotrophic")))) +
    facet_grid(~factor(ecozone, levels = ecozone_longitude_order), scales = "free_x") +
    geom_line(aes(group = lake_id)) +
    labs(x = "Sample size (sequences)", y = "ASV richness") +
    scale_colour_manual(values = palette_trophic, name = "Trophic status") +
    theme_light() %+replace%
    theme(legend.position = "right",
          strip.background = element_rect(fill = "black"),
          axis.text.x = element_text(angle = 330, hjust = 0, vjust = 1)))
#ggsave(filename = "figures/lp2017-2019watercolumn18s_protists_rarefaction_curve_lakes.pdf", plot = rarefaction_curve, device = "pdf", width = 24, height = 3, units = "in")
