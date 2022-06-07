# Assess ASV taxonomic composition and incidence across lakes

# Load libraries
library(tidyverse)
library(ggpubr)
library(patchwork)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import melted ASV data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)

# Import metadata
source("scripts/05_metadata.R")


#### Create taxonomy table ####
taxonomy <- asv_melt %>%
  distinct(asv_code, kingdom, supergroup, division, class, order, family, genus, species)


#### Assess taxon abundances ####
# Calculate taxon contributions to landscape abundance (percentages of sequences)
asv_melt %>%
  group_by(division) %>%
  summarize(nseqs = sum(nseqs)) %>%
  mutate(pct_seqs = nseqs/sum(nseqs) * 100) %>%
  arrange(-pct_seqs)

# Calculate ASV taxonomic composition (ASV richness by group)
taxonomy %>%
  group_by(division) %>%
  dplyr::count(name = "n_asvs") %>%
  ungroup() %>%
  mutate(pct_asvs = n_asvs/sum(n_asvs) * 100) %>%
  arrange(-pct_asvs)

# Calculate mean ASV relative abundance
asvs_mean_relseqs <- asv_melt %>%
  group_by(asv_code) %>%
  summarize(mean_relseqs = mean(relseqs))

# Calculate maximum ASV relative abundance
asvs_max_relseqs <- asv_melt %>%
  group_by(asv_code) %>%
  summarize(max_relseqs = max(relseqs))


#### Assess taxon occurrences ####
(total_lakes <- length(unique(asv_melt$lake_id)))  # n lakes

# Calculate taxon ASV occurrence
# Bin % lakes categories
asv_occurrence <- asv_melt %>%
  group_by(asv_code) %>%
  dplyr::count(name = "n_lakes") %>%
  mutate(pct_lakes = n_lakes/total_lakes * 100) %>%
  mutate(pct_lakes_bin = round(pct_lakes))

# Calculate sequence counts per ASVs
asvs_nseqs <- asv_melt %>%
  group_by(asv_code) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup()

# Assess the relationship between occurrence and mean abundance (for each ASV)
asv_occurrence_relseqs <- asv_occurrence %>%
  left_join(asvs_mean_relseqs, by = "asv_code")

cor(asv_occurrence_relseqs$pct_lakes, asv_occurrence_relseqs$mean_relseqs)

asv_occurrence_relseqs %>%
  left_join(taxonomy, by = "asv_code") %>%
  ggplot(aes(x = pct_lakes, y = mean_relseqs,
             colour = factor(division, levels = c("Ochrophyta", "Opalozoa", "Pseudofungi", "Sagenista", "Stramenopiles_X",
                                                            #"Metazoa", "Fungi"
                                                            "Cryptophyta", "Katablepharidophyta", "Centroheliozoa", "Haptophyta", "Telonemia",
                                                            "Dinoflagellata", "Ciliophora", "Apicomplexa", "Perkinsea", "Protalveolata_X",
                                                            "Chlorophyta", "Streptophyta", "Rhodophyta", "Glaucophyta",
                                                            "Cercozoa", "Foraminifera",
                                                            "Choanoflagellida", "Mesomycetozoa",  # Opisthokonts
                                                            "Discoba", "Metamonada", "Malawimonadidae",
                                                            "Lobosa", "Conosa", "Breviatea",
                                                            "Hilomonadea", "Apusomonadidae")))) +
  geom_point() +
  scale_colour_manual(values = palette_division, na.value = "black") +
  theme_classic() %+replace%
  theme(legend.position = "none")

# Assess the relationship between occurrence and maximum abundance (for each ASV)
asv_occurrence_relseqs <- asv_occurrence_relseqs %>%
  left_join(asvs_max_relseqs, by = "asv_code")

cor(asv_occurrence_relseqs$pct_lakes, asv_occurrence_relseqs$max_relseqs)

# Plot histograms of ASV occurrences (% lakes)
(asvs_lakes_histogram_divisions_part1 <- asv_occurrence %>%
    left_join(taxonomy, by = "asv_code") %>%
    ggplot() +
    geom_histogram(aes(x = pct_lakes, y = ..count..,
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
                   binwidth = 1) +
    xlab("Incidence (% lakes)") +
    ylab("Number of ASVs") +
    scale_fill_manual(values = palette_division, na.value = "black") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(labels = scales::comma, expand = c(0,0)) +
    theme_classic() %+replace%
    theme(legend.position = "none",
          axis.line = element_blank()))

(asvs_lakes_histogram_divisions_part2 <- asv_occurrence %>%
    left_join(taxonomy, by = "asv_code") %>%
    filter(pct_lakes >= 30) %>%
    ggplot() +
    geom_histogram(aes(x = pct_lakes, y = ..count..,
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
                   binwidth = 1) +
    xlab("Occurrence (number of lakes)") +
    ylab("Number of ASVs") +
    scale_fill_manual(values = palette_division, na.value = "black") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(labels = scales::comma, expand = c(0,0)) +
    theme_classic() %+replace%
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.title = element_blank()))

(asvs_lakes_histogram_divisions_all <- asvs_lakes_histogram_divisions_part1 +
    inset_element(asvs_lakes_histogram_divisions_part2,
                  left = 0.2, bottom = 0.2, right = 1, top = 0.6))
#ggsave("figures/lp2017-2019watercolumn18s_occurrence_histogram.pdf", asvs_lakes_histogram_divisions_all, "pdf", width = 6, height = 5, units = "in")

# Plot occurrences of n sequences
(asvs_lakes_histogram_divisions_nseqs <- asv_occurrence %>%
    left_join(asvs_nseqs, by = "asv_code") %>%
    left_join(taxonomy, by = "asv_code") %>%
    group_by(pct_lakes, division) %>%
    summarize(nseqs = sum(nseqs)) %>%
    ggplot() +
    geom_bar(aes(x = pct_lakes, y = nseqs,
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
             stat = "identity") +
    xlab("Incidence (% lakes)") +
    ylab("Number of sequences") +
    scale_fill_manual(values = palette_division, na.value = "black") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() %+replace%
    theme(legend.position = "none",
          axis.line = element_blank()))
#ggsave("figures/lp2017-2019watercolumn18s_occurrence_nseqs_histogram.pdf", asvs_lakes_histogram_divisions_nseqs, "pdf", width = 5, height = 5, units = "in")


#### Taxon accumulations ####
# Rank accumulation curve
# X-axis: Individual ASVs ranked from most to least abundant
# Y-axis: Cumulative percent composition of dataset (in n sequences)
# Overlay ASV occurrence (across lakes) information
percent_comp <- asv_melt %>%
  group_by(asv_code) %>%
  summarize(sum_nseqs = sum(nseqs)) %>%
  mutate(percent_seqs = sum_nseqs/sum(asv_melt$nseqs)*100) %>%
  arrange(desc(percent_seqs), asv_code) %>%
  mutate(percent_accumulation = cumsum(percent_seqs))

(accumulation_plot <- percent_comp %>%
    left_join(asv_occurrence, by = "asv_code") %>%
    mutate(asv_count = 1:n()) %>%
    ggplot(aes(x = asv_count, y = percent_accumulation, colour = pct_lakes)) +
    geom_point(size = 1) +
    xlab("n ASVs") +
    ylab("Accumulated sequence composition (%)") +
    labs(colour = "Occurrence\n(% lakes)") +
    scale_x_continuous(breaks = seq(0, length(unique(percent_comp$asv_code)), 1000)) +
    scale_y_continuous(breaks = seq(0, 100, 5)) +
    scale_colour_gradient2(low = "#24D3FF", mid = "grey", high = "red", midpoint = 50) +
    theme_classic() %+replace%
    theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1)))

# Accumulated ASV count at accumulated sequence compositions
percent_comp_asvs <- percent_comp %>%
  mutate(percent_accumulation_ceiling = ceiling(percent_accumulation)) %>%
  group_by(percent_accumulation_ceiling) %>%
  tally(name = "nasvs") %>%
  arrange(percent_accumulation_ceiling) %>%
  mutate(nasvs_accumulated = cumsum(nasvs))

(asvs_accumulation_plot <- percent_comp_asvs %>%
    ggplot(aes(x = percent_accumulation_ceiling, y = nasvs_accumulated)) +
    geom_point() +
    scale_y_log10() +
    xlab("Accumulated sequence composition (%)") +
    ylab("Accumulated ASV count") +
    theme_classic())

# Combine plots for dataset accumulation
# Are plots #2 and #3 redundant?
(accumulation_plots <- ggarrange(asvs_lakes_histogram, accumulation_plot, asvs_accumulation_plot))
#ggsave("figures/lp2017-2019watercolumn18s_accumulation_plots.pdf", accumulation_plots, "pdf", width = 12, height = 12, units = "in")
