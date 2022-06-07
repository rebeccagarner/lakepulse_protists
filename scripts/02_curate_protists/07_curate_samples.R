# Curate samples

# Load libraries
library(tidyverse)
library(seqinr)


#### Import and format data ####
# Import melted sequence table
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_382lakes.tsv", col_names = TRUE)
length(unique(asv_melt$asv_code))  # Number of ASVs
length(unique(asv_melt$lake_id))  # Number of lakes


#### Remove saline and ELA lake samples ####
# Pull out lake IDs
samples_lakes <- asv_melt %>%
  distinct(lake_id)

# Remove experimental lakes in the ELA
# L223 (06-306) and L302 (06-305) were acidified
# L226 (06-304) L227 (06-307) were nutrient-enriched
# L224 (06-311) and L373 (06-308) are reference lakes
experimental_lakes <- c("06-304", "06-305", "06-306", "06-307")

# Import list of saline lakes (note that not all saline lakes are present in full dataset)
saline_lakes <- read_tsv("output/env/lakes_saline.tsv", col_names = "lake_id") %>%
  pull(lake_id)

# Create lake list
lakes <- samples_lakes %>%
  filter(!lake_id %in% experimental_lakes & !lake_id %in% saline_lakes) %>%
  pull(lake_id)

# Filter melted ASV table so it only contains samples of interest
asv_melt <- asv_melt %>%
  filter(lake_id %in% lakes)
length(unique(asv_melt$lake_id))


#### Remove samples with low sequence counts ####
# Plot histogram of sequence counts by sample
asv_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ggplot() +
  geom_histogram(aes(x = nseqs, y = ..count..), binwidth = 1000) +
  scale_x_continuous() +
  theme_bw()

# Identify samples with fewer than 10,000 sequences
lakes_rm <- asv_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs = sum(nseqs)) %>%
  filter(nseqs < 10000) %>%
  pull(lake_id)

# Filter out low-sequence samples
asv_melt <- asv_melt %>%
  filter(!lake_id %in% lakes_rm)
length(unique(asv_melt$lake_id))


#### Write curated data to files ####
length(unique(asv_melt$lake_id))  # Number of lakes
length(unique(asv_melt$asv_code))  # Number of ASVs
sum(asv_melt$nseqs)  # Total sequences

# Write lake list to file
# asv_melt %>%
#   distinct(lake_id) %>%
#   arrange(lake_id) %>%
#   write_tsv("output/env/lakes.tsv", col_names = FALSE)

# Write new melted ASV table
# asv_melt %>%
#   write_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv", col_names = TRUE)

# Write new fasta file
asvs_seqs <- asv_melt %>%
  distinct(asv_code, sequence) %>%
  arrange(asv_code)

# write.fasta(sequences = as.list(asvs_seqs$sequence),
#             names = asvs_seqs$asv_code,
#             file.out = "output/fasta/lp2017-2019watercolumn18s_protists_366lakes.fasta")
