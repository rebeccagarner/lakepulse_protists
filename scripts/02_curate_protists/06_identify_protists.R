# Identify ASV data assigned to protists

# Load libraries
library(tidyverse)
library(seqinr)
library(ggpubr)
library(ape)
library(ggtree)

# Load palettes
source("scripts/00_palettes.R")


#### Load and format data ####
# Import melted ASV sequence count and taxonomy data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_rmprimers.tsv", col_names = TRUE)
length(unique(asv_melt$asv_code))

# Format taxonomy information
taxonomy <- asv_melt %>%
  distinct(asv_code, sequence, kingdom, supergroup, division, class, order, family, genus, species)


#### Visualize tree taxonomy ####
# # Import tree
# tree_all <- read.tree("output/trees/lp2017-2019watercolumn18s_rmseqs_trimprimers_gtr.tree")
# length(tree_all$tip.label) == length(unique(asv_melt$asv_code))  # Should evaluate to TRUE
# 
# # Join tree and taxonomy
# ggtree_all <- ggtree(tree_all, layout = "rectangular", size = 0.05)
# 
# ggtree_all$data <- ggtree_all$data %>%
#   rename(asv_code = label) %>%
#   left_join(taxonomy, by = "asv_code")
# 
# # Visualize tree with taxonomy labels
# (tree_taxonomy <- ggtree_all +
#     geom_tippoint(data = ggtree_all$data,
#                   aes(colour = division), size = 0.5) +
#     geom_tiplab(data = ggtree_all$data,
#                 aes(label = paste(supergroup, division, class, order, family, genus, sep = "|")), size = 0.2) +
#     scale_colour_manual(values = palette_division, na.value = "black") +
#     theme(legend.position = "none"))
# #ggsave("figures/lp2017-2019watercolumn18s_rmseqs_trimprimers_gtr.pdf", tree_taxonomy, "pdf", width = 12, height = 100, units = "in", limitsize = FALSE)
# 
# (tree_taxonomy_macrobes <- ggtree_all +
#     geom_tippoint(data = ggtree_all$data %>%
#                     filter(division %in% c("Metazoa", "Fungi") | class == "Embryophyceae"),
#                   aes(colour = division), size = 0.5) +
#     geom_tiplab(data = ggtree_all$data,
#                 aes(label = paste(supergroup, division, class, order, family, genus, sep = "|")), size = 0.2) +
#     scale_colour_manual(values = palette_division, na.value = "black") +
#     theme(legend.position = "none"))
# #ggsave("figures/lp2017-2019watercolumn18s_rmseqs_trimprimers_gtr_macrobes.pdf", tree_taxonomy_macrobes, "pdf", width = 12, height = 100, units = "in", limitsize = FALSE)


#### Extract protists ####
# Separate opisthokonts (metazoa and fungi), land plants, and other eukaryotes (protists)
# Metazoa ----
metazoa_melt <- asv_melt %>%
  filter(division == "Metazoa")
metazoa_asvs <- metazoa_melt %>%
  distinct(asv_code) %>%
  pull()

metazoa_seqs <- metazoa_melt %>%
  distinct(asv_code, sequence)
# write.fasta(sequences = as.list(metazoa_seqs$sequence),
#             names = metazoa_seqs$asv_code,
#             file.out = "output/metazoa/lp2017-2019watercolumn18s_metazoa.fasta")

# taxonomy %>%
#   filter(asv_code %in% metazoa_asvs) %>%
#   write_tsv("output/metazoa/lp2017-2019watercolumn18s_metazoa_taxonomy.tsv", col_names = TRUE)

# metazoa_melt %>%
#   select(lake_id, asv_code, nseqs) %>%
#   pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0) %>%
#   write_tsv("output/metazoa/lp2017-2019watercolumn18s_metazoa_sequencetable.tsv", col_names = TRUE)

# Fungi ----
fungi_melt <- asv_melt %>%
  filter(division == "Fungi")
fungi_asvs <- fungi_melt %>%
  distinct(asv_code) %>%
  pull()

fungi_seqs <- fungi_melt %>%
  distinct(asv_code, sequence)
# write.fasta(sequences = as.list(fungi_seqs$sequence),
#             names = fungi_seqs$asv_code,
#             file.out = "output/fungi/lp2017-2019watercolumn18s_fungi.fasta")

# taxonomy %>%
#   filter(asv_code %in% fungi_asvs) %>%
#   write_tsv("output/fungi/lp2017-2019watercolumn18s_fungi_taxonomy.tsv", col_names = TRUE)

# fungi_melt %>%
#   select(lake_id, asv_code, nseqs) %>%
#   pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0) %>%
#   write_tsv("output/fungi/lp2017-2019watercolumn18s_fungi_sequencetable.tsv", col_names = TRUE)

# Land plants ----
plants_melt <- asv_melt %>%
  filter(class == "Embryophyceae")
plant_asvs <- plants_melt %>%
  distinct(asv_code) %>%
  pull()

# Consider remaining ASVs as "protists"
protists_melt <- asv_melt %>%
  filter(!asv_code %in% metazoa_asvs & !asv_code %in% fungi_asvs & !asv_code %in% plant_asvs)

# Calculate relative abundance of protist sequences
lakes_nseqs <- protists_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

protists_melt <- protists_melt %>%
  left_join(lakes_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  select(-nseqs_total)
sum(protists_melt$relseqs)  # Should evaluate to the number of lakes


#### Write melted protists dataset to file ####
length(unique(asv_melt$lake_id))  # Number of lakes
# protists_melt %>%
#   write_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_382lakes.tsv", col_names = TRUE)


#### Tax group ASV/sequence compositions (pie charts) ####
# Pie chart n ASVs assigned to metazoa/fungi/plants/other
(taxgroup_asvs_piechart <- asv_melt %>%
   select(-lake_id, -nseqs) %>%
   distinct(asv_code, .keep_all = TRUE) %>%
   mutate(tax_group = case_when(division == "Metazoa" ~ "metazoa",
                                division == "Fungi" ~ "fungi",
                                class == "Embryophyceae" ~ "plants",
                                TRUE ~ "protists")) %>%
   group_by(tax_group) %>%
   tally(name = "nasvs") %>%
   ggplot(aes(x = "", y = nasvs, fill = tax_group)) +
   geom_bar(stat = "identity") +
   coord_polar("y") +
   geom_text(aes(label = str_c(nasvs, " (", round(nasvs/(length(unique(asv_melt$asv_code)))*100, 1), "%)"), x = 1),
             position = position_stack(vjust = 0.5),
             fontface = "bold", size = 4) +
   scale_fill_manual(values = palette_taxgroup) +
   theme_void())

# Pie chart n sequences assigned to metazoa/fungi/plants/other
(taxgroup_seqs_piechart <- asv_melt %>%
    select(-lake_id) %>%
    mutate(tax_group = case_when(division == "Metazoa" ~ "metazoa",
                                 division == "Fungi" ~ "fungi",
                                 class == "Embryophyceae" ~ "plants",
                                 TRUE ~ "protists")) %>%
    group_by(tax_group) %>%
    summarize(sum_nseqs = sum(nseqs)) %>%
    ggplot(aes(x = "", y = sum_nseqs, fill = tax_group)) +
    geom_bar(stat = "identity") +
    coord_polar("y") +
    geom_text(aes(label = str_c(sum_nseqs, " (", round(sum_nseqs/sum(asv_melt$nseqs)*100, 1), "%)"), x = 1),
              position = position_stack(vjust = 0.5),
              fontface = "bold", size = 4) +
    scale_fill_manual(values = palette_taxgroup) +
    theme_void())

# Combine pie charts n ASVs and n sequences assigned to metazoa/fungi/plants/other
(taxgroup_piecharts <- ggarrange(taxgroup_asvs_piechart, taxgroup_seqs_piechart,
                                 labels = c("ASVs", "Sequences"),
                                 common.legend = TRUE, legend = "right"))
#ggsave("figures/lp2017-2019watercolumn18s_taxgroup_piecharts.pdf", taxgroup_piecharts, "pdf", width = 10, height = 5, units = "in")
