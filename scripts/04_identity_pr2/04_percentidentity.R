# Assess sequence identities between ASVs and PR2 references

# Load libraries
library(tidyverse)
library(ggpubr)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import melted ASV data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)

# Import BLAST output:
# - qseqid    Query Seq-id
# - sseqid    Subject Seq-id
# - pident    Percentage of identical matches
# - length    Alignment length
# - mismatch  Number of mismatches
# - gapopen   Number of gap openings
# - qlen      Query sequence length
# - qstart    Start of alignment in query
# - qend      End of alignment in query
# - slen      Subject sequence length
# - sstart    Start of alignment in subject
# - send      End of alignment in subject
# - evalue    Expect value
# - bitscore  Bit score
blast <- read_tsv("output/pr2/lp2017-2019watercolumn18s_protists_366lakes_blast_pr2v7.tsv",
                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                "gapopen", "qlen", "qstart", "qend", "slen",
                                "sstart", "send", "evalue", "bitscore"))

# Format taxonomy in BLAST output
blast_taxonomy <- blast %>%
  mutate(sseqid = str_remove(sseqid, ";$")) %>%
  separate(sseqid, into = c("domain", "supergroup", "division", "class", "family", "order", "genus", "species"), sep = ";")


#### Isolate globally aligned queries ####
# Extract top hits
blast_tophits <- blast_taxonomy %>%
  group_by(qseqid) %>%
  slice_max(order_by = bitscore) %>%
  ungroup()
length(unique(blast_taxonomy$qseqid)) == nrow(blast_tophits)  # Should evaluate to TRUE

# Retain BLAST results for ASVs aligned over >= 220 nucleotides
blast_global <- blast_tophits %>%
  filter(length >= 220)
length(unique(blast_global$qseqid))  # Number of globally aligned ASVs


#### Bin sequence identities ####
blast_binned <- blast_global %>%
  mutate(identity_bin = round(pident))


#### Visualize sequence identities distribution ####
# Extend taxonomic division colour palette to include other taxa in PR2
palette_division_other <- c("black", "black", "black", "black", "black",
                            "black", "black", "black", "black")
names(palette_division_other) <- c("Radiolaria", "Picozoa", "Alveolata_X", "Opisthokonta_X", "Amoebozoa_X",
                                   "Eukaryota_XX", "Mantamonadidea", "Hacrobia_X", "Prasinodermophyta")

palette_division <- c(palette_division, palette_division_other)

# Plot n ASVs distributed across sequence identities
(blast_nasvs_histogram <- blast_binned %>%
    group_by(identity_bin, division) %>%
    dplyr::count(name = "n_asvs") %>%
    ggplot() +
    geom_bar(aes(x = identity_bin, y = n_asvs,
                       fill = factor(division, levels = c("Ochrophyta", "Opalozoa", "Pseudofungi", "Sagenista", "Stramenopiles_X",
                                                          "Metazoa", "Fungi", "Choanoflagellida", "Mesomycetozoa",
                                                          "Cryptophyta", "Katablepharidophyta", "Centroheliozoa", "Haptophyta", "Telonemia",
                                                          "Dinoflagellata", "Ciliophora", "Apicomplexa", "Perkinsea", "Protalveolata_X",
                                                          "Chlorophyta", "Streptophyta", "Rhodophyta", "Glaucophyta",
                                                          "Cercozoa", "Foraminifera",
                                                          "Discoba", "Metamonada", "Malawimonadidae",
                                                          "Lobosa", "Conosa", "Breviatea",
                                                          "Hilomonadea", "Apusomonadidae",
                                                          
                                                          "Radiolaria", "Picozoa", "Alveolata_X", "Opisthokonta_X", "Amoebozoa_X", "Eukaryota_XX", "Mantamonadidea", "Hacrobia_X", "Prasinodermophyta"))),
             stat = "identity", width = 0.9) +
    scale_fill_manual(values = palette_division, na.value = "lightgrey", name = "PR2 taxonomy") +
    scale_x_reverse(breaks = seq(0, 100, by = 5)) +
    labs(x = "Sequence identity (%)",
         y = "Number of ASVs") +
    theme_classic() %+replace%
    theme(legend.position = "right",
          panel.background = element_blank(), panel.grid = element_blank()))

# Calculate sequence counts per sequence identity category
asvs_nseqs <- asv_melt %>%
  group_by(asv_code) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup()

# Plot n sequences distributed across sequence identities
(blast_nseqs_histogram <- blast_binned %>%
    left_join(asvs_nseqs, by = c("qseqid" = "asv_code")) %>%
    group_by(identity_bin, division) %>%
    summarize(nseqs = sum(nseqs)) %>%
    ggplot() +
    geom_bar(aes(x = identity_bin, y = nseqs,
                 fill = factor(division, levels = c("Ochrophyta", "Opalozoa", "Pseudofungi", "Sagenista", "Stramenopiles_X",
                                                    "Metazoa", "Fungi", "Choanoflagellida", "Mesomycetozoa",
                                                    "Cryptophyta", "Katablepharidophyta", "Centroheliozoa", "Haptophyta", "Telonemia",
                                                    "Dinoflagellata", "Ciliophora", "Apicomplexa", "Perkinsea", "Protalveolata_X",
                                                    "Chlorophyta", "Streptophyta", "Rhodophyta", "Glaucophyta",
                                                    "Cercozoa", "Foraminifera",
                                                    "Discoba", "Metamonada", "Malawimonadidae",
                                                    "Lobosa", "Conosa", "Breviatea",
                                                    "Hilomonadea", "Apusomonadidae",
                                                    
                                                    "Radiolaria", "Picozoa", "Alveolata_X", "Opisthokonta_X", "Amoebozoa_X", "Eukaryota_XX", "Mantamonadidea", "Hacrobia_X", "Prasinodermophyta"))),
             stat = "identity", width = 0.9) +
    scale_fill_manual(values = palette_division, na.value = "lightgrey", name = "PR2 taxonomy") +
    scale_x_reverse(breaks = seq(0, 100, by = 5)) +
    labs(x = "Sequence identity (%)",
         y = "Number of sequences") +
    theme_classic() %+replace%
    theme(legend.position = "right",
          panel.background = element_blank(), panel.grid = element_blank()))

# Combine sequence identity plots
(blast_histogram_all <- ggarrange(blast_nasvs_histogram, blast_nseqs_histogram, nrow = 1, common.legend = TRUE, legend = "right"))
#ggsave("figures/lp2017-2019watercolumn18s_pr2_identity_histograms.pdf", blast_histogram_all, "pdf", width = 12, height = 5, units = "in")

# Calculate ASV occurrence
occurrence <- asv_melt %>%
  group_by(asv_code) %>%
  count(name = "n_lakes")

blast_binned %>%
  left_join(asvs_nseqs, by = c("qseqid" = "asv_code")) %>%
  left_join(occurrence, by = c("qseqid" = "asv_code")) %>%
  #group_by(identity_bin, division) %>%
  #summarize(nseqs = sum(nseqs)) %>%
  ggplot() +
  geom_point(aes(x = n_lakes, y = nseqs, colour = pident)) +
  scale_y_log10() +
  scale_color_gradient2(low = "red", mid = "white", high = "#2A7CAB", midpoint = 85) +
  labs(x = "Incidence (number of lakes)",
       y = "Number of sequences",
       colour = "Identity with reference (%)") +
  #theme_classic() %+replace%
  theme(legend.position = "right",
        panel.background = element_blank(), panel.grid = element_blank())

# Assess ASV frequencies at various sequence identity thresholds
identity_nasvs <- tibble(identity = NA,
                         n_asvs = NA)

identity_min <- floor(min(blast_global$pident))
identity_max <- ceiling(max(blast_global$pident))
identity <- identity_min

while (identity <= identity_max) {
  blast_global_tmp <- blast_global %>%
    filter(pident >= identity)
  n_asvs <- nrow(blast_global_tmp)
  identity_tmp <- tibble(identity, n_asvs)
  identity_nasvs <- bind_rows(identity_nasvs, identity_tmp)
  identity <- identity + 1
  print(identity)
}

# Calculate n sequences at 100% sequence identity
blast_binned %>%
  left_join(asvs_nseqs, by = c("qseqid" = "asv_code")) %>%
  filter(pident == 100) %>%
  summarize(nseqs = sum(nseqs))


#### Assess taxonomic assignments of novel genotypes ####
# Create a taxonomy table
taxonomy <- asv_melt %>%
  distinct(asv_code, supergroup, division, class, family, order, genus, species)

# Count ASVs with low sequence similarity assigned to various taxonomic groups
blast_binned %>%
  left_join(asvs_nseqs, by = c("qseqid" = "asv_code")) %>%
  select(-domain, -supergroup, -division, -class, -family, -order, -genus, -species) %>%
  filter(pident < 90) %>%
  left_join(taxonomy, by = c("qseqid" = "asv_code")) %>%
  group_by(supergroup) %>%
  dplyr::count() %>%
  arrange(-n)
