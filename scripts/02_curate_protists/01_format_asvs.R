# Join ASV (sequence and taxonomy) data

# Load libraries
library(tidyverse)
library(seqinr)


#### Load and format data ####
# Load DADA2 sequence table
load("output/dada2/lp2017-2019watercolumn18s_seqtab_nochim.rda")

# Load taxonomy table of PR2 classifications
load("output/dada2/lp2017-2019watercolumn18s_taxonomy_pr2v412.rda")

taxonomy <- taxonomy %>%
  as_tibble(rownames = "sequence")


#### Assign ASV codes ####
# Number of unique ASVs
(nasvs <- n_distinct(taxonomy$sequence))

# Assign ASV codes to unique sequences
taxonomy <- taxonomy %>%
  mutate(asv_code = paste0("ASV", str_pad(string = 1:nasvs, width = nchar(nasvs), side = "left", pad = "0")))


#### Write sequences to fasta file ####
# write.fasta(sequences = as.list(taxonomy$sequence),
#             names = taxonomy$asv_code,
#             file.out = "output/fasta/lp2017-2019watercolumn18s_all.fasta")


#### Combine sequence counts and taxonomy ####
# Convert sequence table to long format and join taxonomy
asv_melt <- seqtab.nochim %>%
  as_tibble(rownames = "sample_id") %>%
  mutate(lake_id = str_remove(string = sample_id, pattern = "LP-watercolumn-18S_")) %>%
  select(-sample_id) %>%
  pivot_longer(!lake_id, names_to = "sequence", values_to = "nseqs") %>%
  filter(nseqs > 0) %>%
  left_join(taxonomy, by = "sequence")

# Write melted sequence data to file
# asv_melt %>%
#   write_tsv("output/melted/lp2017-2019watercolumn18s_melt_all.tsv", col_names = TRUE)
