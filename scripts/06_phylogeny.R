# Plot phylogenetic tree for protist ASVs

# Load libraries
library(tidyverse)
library(seqinr)
library(ggpubr)
library(ape)
library(ggtree)

# Load palettes
source("scripts/00_palettes.R")


#### Load and format ASV data ####
# Import melted protists data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)
length(unique(asv_melt$lake_id))  # Number of lakes

# Import rooted phylogenetic tree
tree_protists <- read.tree("output/trees/lp2017-2019watercolumn18s_protists_366lakes_sina_gtr_rooted.tree")
length(tree_protists$tip.label) == length(unique(asv_melt$asv_code))  # Should evaluate to TRUE

# Format taxonomy information
taxonomy <- asv_melt %>%
  distinct(asv_code, sequence, kingdom, supergroup, division, class, order, family, genus, species)


#### Plot circular tree with branches coloured by taxonomy ####
ggtree_protists <- ggtree(tree_protists, layout = "circular", size = 0.1)  # Size is branch line width
ggtree_protists$data <- ggtree_protists$data %>%
  rename(asv_code = label) %>%
  left_join(taxonomy, by = "asv_code") %>%
  mutate(division = case_when(is.na(division) ~ "unclassified",
                              TRUE ~ division))

# Format tree
tree_data_tips <- ggtree_protists$data %>%
  filter(isTip)

divisions_list <- ggtree_protists$data %>%
  filter(isTip) %>%
  select(division, asv_code) %>%
  group_split(division)

divisions_list %>%
  purrr::map(~pull(.,division)) %>%  # Pull out division variable
  purrr::map(~as.character(.)) %>%  # Convert factor to character
  purrr::map(~unique(.)) -> names(divisions_list)

divisions_list <- divisions_list %>%
  map(~select(., -division)) %>%
  map(~pull(., asv_code))

tree_protists <- groupOTU(tree_protists, divisions_list, group_name = "division")

(tree_all <- ggtree(tree_protists,  aes(colour = division), layout = "circular", size = 0.1) +
    #geom_tippoint(data = ggtree_protists$data, aes(colour = division), size = 0.5) +
    labs(colour = "division") +
    scale_colour_manual(values = palette_division, na.value = "black") +
    theme(legend.position = "none"))
#ggsave("figures/lp2017-2019watercolumn18s_protists_366lakes_sina_gtr_rooted.pdf", tree_all, "pdf", width = 30, height = 30, units = "in", limitsize = FALSE)
