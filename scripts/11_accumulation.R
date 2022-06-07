# Plot ASV accumulation curve

# Load libraries
library(tidyverse)
library(vegan)


#### Import and format data ####
# Import melted protists data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)


#### Create ASV tables ####
# ASV table of sequence counts
asvtable_nseqs <- asv_melt %>%
  select(lake_id, asv_code, nseqs) %>%
  spread(asv_code, nseqs, fill = 0) %>%
  column_to_rownames("lake_id")


#### Species accumulation curve ####
# Compute ASV accumulation
(accumulation <- specaccum(asvtable_nseqs, method = "random", permutations = 100))

# Summarize ASV accumulation
accumulation_summary <- tibble(n_lakes = accumulation$sites,
                               n_asvs = accumulation$richness,
                               sd = accumulation$sd)

# Write ASV accumulation summary to file
# accumulation_summary %>%
#   write_tsv("output/diversity/lp2017-2019watercolumn18s_accumulation.tsv", col_names = TRUE)

# Import accumulation data
accumulation_summary <- read_tsv("output/diversity/lp2017-2019watercolumn18s_accumulation.tsv", col_names = TRUE)

# Plot accumulation curve
(accumulation_curve <- accumulation_summary %>%
    ggplot(aes(x = n_lakes, y = n_asvs)) +
    geom_point(size = 1) +
    geom_linerange(aes(ymin = n_asvs - sd, ymax = n_asvs + sd),
                    alpha = 0.2) +
    labs(x = "Number of lakes sampled",
         y = "ASV richness") +
    scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
    scale_y_continuous(labels = scales::comma, expand = c(0,0)) +
    theme_light() %+replace%
    theme(legend.position = "none",
          panel.border = element_blank()))
#ggsave(filename = "figures/lp2017-2019watercolumn18s_accumulation_curve.pdf", plot = accumulation_curve, device = "pdf", width = 6, height = 5, units = "in")
