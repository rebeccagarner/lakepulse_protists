# Assess PR2 18S rRNA gene V7 region fragment lengths

# Load libraries
library(tidyverse)
library(seqinr)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import PR2 V7 region fragments
pr2_v7 <- read.fasta("output/pr2/pr2_version_4.13.0_18S_v7_dada2.fasta", "DNA", as.string = TRUE, forceDNAtolower = FALSE)
pr2_v7 <- tibble(header = names(pr2_v7),
                 sequence = unlist(getSequence(pr2_v7, as.string = TRUE)))

# Format taxonomy in PR2 sequence headers
pr2_v7_taxonomy <- pr2_v7 %>%
  mutate(header = str_remove(header, ";$")) %>%
  separate(header, into = c("domain", "supergroup", "division", "class", "family", "order", "genus", "species"), sep = ";")


#### Assess V7 region fragment lengths ####
# Calculate sequence lengths
pr2_v7_seqlength <- pr2_v7_taxonomy %>%
  mutate(seq_length = nchar(sequence))

mean(pr2_v7_seqlength$seq_length)  # Mean sequence length

# Assess sequence length for putative protist sequences
pr2_v7_seqlength %>%
  filter(division != "Metazoa" & division != "Fungi" & class != "Embryophyceae") %>%
  summarize(mean_length = mean(seq_length))

# Extend taxonomic division colour palette to include other taxa in PR2
palette_division_other <- c("black", "black", "black", "black", "black",
                            "black", "black", "black", "black")
names(palette_division_other) <- c("Radiolaria", "Picozoa", "Alveolata_X", "Opisthokonta_X", "Amoebozoa_X",
                                   "Eukaryota_XX", "Mantamonadidea", "Hacrobia_X", "Prasinodermophyta")

palette_division <- c(palette_division, palette_division_other)

# Plot sequence length distribution
pr2_v7_seqlength %>%
  ggplot() +
  geom_histogram(aes(x = seq_length, y = ..count..,
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
                 binwidth = 10) +
  scale_fill_manual(values = palette_division, na.value = "lightgrey", name = "Division") +
  scale_x_continuous(breaks = seq(0, 1800, by = 100)) +
  labs(x = "Sequence length (nt)",
       y = "Number of reference sequences") +
  theme_classic() %+replace%
  theme(legend.position = "right")
