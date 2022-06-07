# Remove remaining primer sequences

# Load libraries
library(tidyverse)
library(seqinr)


#### Import and format data ####
# Import melted ASV data
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_rmpos.tsv", col_names = TRUE)
length(unique(asv_melt$asv_code)) == length(unique(asv_melt$sequence))  # Does the number of ASVs equate to the number of sequences?


#### Calculate ASV sequence counts ####
asv_nseqs <- asv_melt %>%
  group_by(asv_code) %>%
  summarize(nseqs = sum(nseqs))


#### Import new alignment ####
alignment <- read.fasta("output/fasta/lp2017-2019watercolumn18s_rmpos_mafft.fasta", "DNA", as.string = TRUE, forceDNAtolower = FALSE)
alignment <- tibble(asv_code = names(alignment),
                    sequence = unlist(getSequence(alignment, as.string = TRUE)))

asv_start_position <- 3L
asv_end_position <- 923L

# Identify ASVs with retained adapter/primer sequences
asvs_wprimers <- alignment %>%
  mutate(rm_fwdseq = str_sub(sequence, start = 1L, end = asv_start_position - 1L),
         rm_revseq = str_sub(sequence, start = asv_end_position + 1L, end = -1L)) %>%
  mutate(rm_fwdseq = str_remove_all(rm_fwdseq, "-"),
         rm_revseq = str_remove_all(rm_revseq, "-")) %>%
  mutate(trim = case_when(rm_fwdseq != "" | rm_revseq != "" ~ TRUE,
                          TRUE ~ FALSE)) %>%
  filter(trim) %>%
  pull(asv_code)

# Assess n sequences associated with these ASVs
asv_nseqs %>%
  filter(asv_code %in% asvs_wprimers)

# Remove ASVs with adapter/primer sequences from melted ASV table
asv_melt <- asv_melt %>%
  filter(!asv_code %in% asvs_wprimers)
length(unique(asv_melt$asv_code))
length(unique(asv_melt$sequence))


#### Write new files for retained ASVs ####
asv_rmprimers <- asv_melt %>%
  distinct(asv_code, sequence) %>%
  arrange(asv_code)

# Write new fasta file for retained ASVs
# write.fasta(sequences = as.list(asv_rmprimers$sequence),
#             names = asv_rmprimers$asv_code,
#             file.out = "output/fasta/lp2017-2019watercolumn18s_rmprimers.fasta")

# Write new melted ASV table for retained ASVs
# asv_melt %>%
#   write_tsv("output/melted/lp2017-2019watercolumn18s_melt_rmprimers.tsv")
