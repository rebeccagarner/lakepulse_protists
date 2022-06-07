# Plot ASV rarefaction curve (landscape)

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import rarefaction tables
rarefaction_protists <- read_tsv("output/diversity/lp2017-2019watercolumn18s_protists_rarefaction.tsv", col_names = TRUE) %>%
  mutate(taxon = "protists")


#### Plot rarefaction curve ####
(rarefaction_curve <- rarefaction_protists %>%
   ggplot(aes(x = rarefied, y = richness,
              colour = factor(taxon, levels = c("Ochrophyta", "Opalozoa", "Pseudofungi", "Sagenista", "Stramenopiles_X",
                                                #"Metazoa", "Fungi"
                                                "Cryptophyta", "Katablepharidophyta", "Centroheliozoa", "Haptophyta", "Telonemia",
                                                "Dinoflagellata", "Ciliophora", "Apicomplexa", "Perkinsea", "Protalveolata_X",
                                                "Chlorophyta", "Streptophyta", "Rhodophyta", "Glaucophyta",
                                                "Cercozoa", "Foraminifera",
                                                "Choanoflagellida", "Mesomycetozoa",  # Opisthokonts
                                                "Discoba", "Metamonada", "Malawimonadidae",
                                                "Lobosa", "Conosa", "Breviatea",
                                                "Hilomonadea", "Apusomonadidae")))) +
   geom_line() +
   labs(x = "Sample size (sequences)", y = "Number of ASVs") +
   scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
   scale_y_continuous(labels = scales::comma, expand = c(0,0)) +
   scale_colour_manual(values = palette_division, na.value = "black") +
   theme_light() %+replace%
   theme(legend.position = "none",
         panel.border = element_blank()))
#ggsave(filename = "figures/lp2017-2019watercolumn18s_protists_rarefaction_curve.pdf", plot = rarefaction_curve, device = "pdf", width = 6, height = 5, units = "in")
