# Assign trophic functional diversity

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import melted ASV data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)

# Import trophic function annotations
trophic <- read_csv("data/trophic_function/pr2_adl_trophic_functions.csv", col_names = TRUE)


#### Link taxonomy and function ####
# Create taxonomy table
taxonomy <- asv_melt %>%
  distinct(asv_code, kingdom, supergroup, division, class, order, family, genus, species)

# Link taxonomy and trophic function
taxonomy_function <- taxonomy
taxonomy_function$trophic_function <- NA
for (i in 1:nrow(taxonomy_function)) {
  if (!is.na(taxonomy_function$supergroup[i]) & taxonomy_function$supergroup[i] %in% trophic$supergroup) {
    trophic_tmp <- trophic %>%
      filter(!is.na(supergroup) & is.na(division) & is.na(class) & is.na(family) & is.na(order) & is.na(genus) & is.na(species))
    if (taxonomy_function$supergroup[i] %in% trophic_tmp$supergroup) {
      taxonomy_function$trophic_function[i] <- trophic_tmp$trophic_function[which(trophic_tmp$supergroup == taxonomy_function$supergroup[i])]
    }
  }
  if (!is.na(taxonomy_function$division[i]) & taxonomy_function$division[i] %in% trophic$division) {
    trophic_tmp <- trophic %>%
      filter(!is.na(division) & is.na(class) & is.na(order) & is.na(family) & is.na(genus) & is.na(species))
    if (taxonomy_function$division[i] %in% trophic_tmp$division) {
      taxonomy_function$trophic_function[i] <- trophic_tmp$trophic_function[which(trophic_tmp$division == taxonomy_function$division[i])]
    }
  }
  if (!is.na(taxonomy_function$class[i]) & taxonomy_function$class[i] %in% trophic$class) {
    trophic_tmp <- trophic %>%
      filter(!is.na(class) & is.na(order) & is.na(family) & is.na(genus) & is.na(species))
    if (taxonomy_function$class[i] %in% trophic_tmp$class) {
      taxonomy_function$trophic_function[i] <- trophic_tmp$trophic_function[which(trophic_tmp$class == taxonomy_function$class[i])]
    }
  }
  if (!is.na(taxonomy_function$order[i]) & taxonomy_function$order[i] %in% trophic$order) {
    trophic_tmp <- trophic %>%
      filter(!is.na(order) & is.na(family) & is.na(genus) & is.na(species))
    if (taxonomy_function$order[i] %in% trophic_tmp$order) {
      taxonomy_function$trophic_function[i] <- trophic_tmp$trophic_function[which(trophic_tmp$order == taxonomy_function$order[i])]
    }
  }
  if (!is.na(taxonomy_function$family[i]) & taxonomy_function$family[i] %in% trophic$family) {
    trophic_tmp <- trophic %>%
      filter(!is.na(family) & is.na(genus) & is.na(species))
    if (taxonomy_function$family[i] %in% trophic_tmp$family) {
      taxonomy_function$trophic_function[i] <- trophic_tmp$trophic_function[which(trophic_tmp$family == taxonomy_function$family[i])]
    }
  }
  if (!is.na(taxonomy_function$genus[i]) & taxonomy_function$genus[i] %in% trophic$genus) {
    trophic_tmp <- trophic %>%
      filter(!is.na(genus) & is.na(species))
    if (taxonomy_function$genus[i] %in% trophic_tmp$genus) {
      taxonomy_function$trophic_function[i] <- trophic_tmp$trophic_function[which(trophic_tmp$genus == taxonomy_function$genus[i])]
    }
  }
  if (!is.na(taxonomy_function$species[i]) & taxonomy_function$species[i] %in% trophic$species) {
    trophic_tmp <- trophic %>%
      filter(!is.na(species))
    if (taxonomy_function$species[i] %in% trophic_tmp$species) {
      taxonomy_function$trophic_function[i] <- trophic_tmp$trophic_function[which(trophic_tmp$species == taxonomy_function$species[i])]
    }
  }
}

# Recategorize trophic groups
taxonomy_function <- taxonomy_function %>%
  mutate(trophic_group = case_when(trophic_function == "bacterivory" ~ "consumers",
                                   trophic_function == "bacterivory, cytotrophy" ~ "consumers",
                                   trophic_function == "bacterivory, detritivory, cytotrophy, invertebrates" ~ "consumers",
                                   trophic_function == "bacterivory, detritivory, osmotrophy" ~ "consumers",
                                   trophic_function == "bacterivory, mixotrophy" ~ "consumers, mixotrophs",
                                   trophic_function == "bacterivory, parasitism" ~ "consumers, parasites",
                                   trophic_function == "bacterivory, saprotrophy" ~ "consumers",
                                   trophic_function == "commensalism" ~ "parasites",
                                   trophic_function == "commensalism, parasitism" ~ "parasites",
                                   trophic_function == "commensalism, parasitism, bacterivory" ~ "consumers, parasites",
                                   trophic_function == "cytotrophy" ~ "consumers",
                                   trophic_function == "cytotrophy, invertebrates" ~ "consumers",
                                   trophic_function == "cytotrophy, mixotrophy" ~ "consumers, mixotrophs",
                                   trophic_function == "cytotrophy, parasitism, commensalism" ~ "consumers, parasites",
                                   trophic_function == "heterotrophy" ~ "consumers",
                                   trophic_function == "heterotrophy, cytotrophy" ~ "consumers",
                                   trophic_function == "heterotrophy, mixotrophy, bacterivory" ~ "consumers, mixotrophs",
                                   trophic_function == "heterotrophy, parasitism" ~ "consumers, parasites",
                                   trophic_function == "heterotrophy, photoautotrophy" ~ "phototrophs, consumers",
                                   trophic_function == "invertebrates" ~ "consumers",
                                   trophic_function == "mixotrophy" ~ "mixotrophs",
                                   trophic_function == "mixotrophy, bacterivory" ~ "consumers, mixotrophs",
                                   trophic_function == "mixotrophy, bacterivory, cytotrophy" ~ "consumers, mixotrophs",
                                   trophic_function == "mixotrophy, cytotrophy" ~ "consumers, mixotrophs",
                                   trophic_function == "mixotrophy, heterotrophy" ~ "consumers, mixotrophs",
                                   trophic_function == "mycotrophy or phycotrophy" ~ "consumers",
                                   trophic_function == "mycotrophy, bacterivory" ~ "consumers",
                                   trophic_function == "mycotrophy, cytotrophy, bacterivory" ~ "consumers",
                                   trophic_function == "omnivory" ~ "consumers",
                                   trophic_function == "osmotrophy" ~ "consumers",
                                   trophic_function == "parasitism" ~ "parasites",
                                   trophic_function == "parasitism, commensalism" ~ "parasites",
                                   trophic_function == "parasitism, histophagy" ~ "parasites",
                                   trophic_function == "pathogenicity" ~ "parasites",
                                   trophic_function == "photoautotrophy" ~ "phototrophs",
                                   trophic_function == "photoautotrophy, cytotrophy" ~ "phototrophs, consumers",
                                   trophic_function == "photoautotrophy, mixotrophy" ~ "phototrophs, mixotrophs",
                                   trophic_function == "photoautotrophy, mixotrophy, bacterivory" ~ "phototrophs, consumers, mixotrophs",
                                   trophic_function == "photoautotrophy, mixotrophy, cytotrophy" ~ "phototrophs, consumers, mixotrophs",
                                   trophic_function == "photoautotrophy, some mixotrophy" ~ "phototrophs, mixotrophs",
                                   trophic_function == "saprotrophy" ~ "consumers",
                                   trophic_function == "saprotrophy or parasitism" ~ "consumers, parasites",
                                   trophic_function == "symbiosis with phototroph" ~ ""))


#### Assess functional diversity ####
# Join trophic function information with melted ASV data
asv_melt <- asv_melt %>%
  left_join(taxonomy_function, by = c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species", "asv_code"))

# Write functional annotations to melted ASV data file
# asv_melt %>%
#   write_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_function.tsv", col_names = TRUE)

# Write taxonomic and functional annotations for unique ASVs to file
# asv_melt %>%
#   distinct(asv_code, sequence,
#            kingdom, supergroup, division, class, order, family, genus, species,
#            trophic_function, trophic_group) %>%
#   select(asv_code, sequence,
#          kingdom, supergroup, division, class, order, family, genus, species,
#          trophic_function, trophic_group) %>%
#   arrange(asv_code) %>%
#   mutate(trophic_group = str_replace(trophic_group, "consumer", "heterotroph")) %>%
#   write_csv("output/otutables/lp2017-2019watercolumn18s_asv_taxonomy_trophic.csv", col_names = TRUE)

# Calculate n ASVs assigned to each trophic function
taxonomy_function %>%
  group_by(trophic_function) %>%
  dplyr::count(name = "n_asvs") %>%
  arrange(-n_asvs) %>%
  ungroup() %>%
  mutate(pct_asvs = n_asvs/length(unique(taxonomy_function$asv_code)) * 100)

taxonomy_function_rmna <- taxonomy_function %>%
  filter(!is.na(trophic_function))
(nasvs_fct_annotated <- length(unique(taxonomy_function_rmna$asv_code)))  # n ASVs with functional annotations
nasvs_fct_annotated/length(unique(taxonomy_function$asv_code)) * 100  # Percent of ASVs with functional annotations

taxonomy_function_rmna %>%
  group_by(trophic_function) %>%
  dplyr::count(name = "n_asvs") %>%
  arrange(-n_asvs) %>%
  ungroup() %>%
  mutate(pct_asvs = n_asvs/nasvs_fct_annotated * 100)

# Calculate n sequences assigned to each trophic function
asv_melt %>%
  group_by(trophic_function) %>%
  summarize(nseqs = sum(nseqs)) %>%
  arrange(-nseqs) %>%
  mutate(pct_seqs = nseqs/sum(asv_melt$nseqs) * 100)

asv_melt_rmna <- asv_melt %>%
  filter(!is.na(trophic_function))
(nseqs_fct_annotated <- sum(asv_melt_rmna$nseqs))  # n sequences with functional annotations
nseqs_fct_annotated/sum(asv_melt$nseqs) * 100  # Percentage of sequences with functional annotations

asv_melt_rmna %>%
  group_by(trophic_function) %>%
  summarize(nseqs = sum(nseqs)) %>%
  arrange(-nseqs)

# Calculate n ASVs assigned to each trophic group
taxonomy_function_rmna %>%
  group_by(trophic_group) %>%
  dplyr::count(name = "n_asvs") %>%
  arrange(-n_asvs) %>%
  ungroup() %>%
  mutate(pct_asvs = n_asvs/nasvs_fct_annotated * 100)

# Calculate n sequences assigned to each trophic group
asv_melt %>%
  group_by(trophic_group) %>%
  summarize(nseqs = sum(nseqs)) %>%
  arrange(-nseqs) %>%
  mutate(pct_seqs = nseqs/sum(asv_melt$nseqs) * 100)
