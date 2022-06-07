# Trim remaining primer sequence positions

# Load libraries
library(tidyverse)
library(seqinr)


#### Import and format data ####
# Import melted sequence table
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_rmspurious.tsv", col_names = TRUE)
length(unique(asv_melt$asv_code))  # Number of ASVs

# Import alignment
alignment <- read.fasta("output/fasta/lp2017-2019watercolumn18s_rmspurious_sina.fasta", "DNA", as.string = TRUE, forceDNAtolower = FALSE)
alignment <- tibble(asv_code = names(alignment),
                    sequence = unlist(getSequence(alignment, as.string = TRUE)))


#### Calculate ASV sequence counts ####
asv_nseqs <- asv_melt %>%
  group_by(asv_code) %>%
  summarize(nseqs = sum(nseqs))


#### Trim poorly aligned positions ####
# Define start and end positions
asv_start_position <- 398L
asv_end_position <- 2039L

# Identify ASVs aligned outside target region
asvs_outsidepos <- alignment %>%
  mutate(sequence = str_replace_all(sequence, "U", "T")) %>%
  mutate(rm_fwdseq = str_sub(sequence, start = 1L, end = asv_start_position - 3L),
         rm_revseq = str_sub(sequence, start = asv_end_position + 1L, end = -3L)) %>%
  mutate(rm_fwdseq = str_remove_all(rm_fwdseq, "-"),
         rm_revseq = str_remove_all(rm_revseq, "-")) %>%
  mutate(tf_fwdseq = case_when(rm_fwdseq != "" ~ TRUE,
                               TRUE ~ FALSE),
         tf_revseq = case_when(rm_revseq != "" ~ TRUE,
                               TRUE ~ FALSE)) %>%
  filter(tf_fwdseq | tf_revseq) %>%
  pull(asv_code)

# Remove ASVs aligned outside target region
asv_melt <- asv_melt %>%
  filter(!asv_code %in% asvs_outsidepos)
length(unique(asv_melt$asv_code))  # Number of ASVs

# Trim remaining ASVs to target positions
asv_trim <- alignment %>%
  filter(asv_code %in% unique(asv_melt$asv_code)) %>%
  mutate(sequence = str_replace_all(sequence, "U", "T")) %>%
  mutate(trimmed = str_sub(sequence, start = asv_start_position, end = asv_end_position)) %>%
  mutate(trimmed = str_replace_all(trimmed, "-", "")) %>%
  mutate(seq_length = nchar(trimmed))

# Assess fragment length lower and upper limits
asv_trim %>%
  select(asv_code, seq_length) %>%
  ggplot() +
  geom_histogram(aes(x = seq_length, y = ..count..), binwidth = 1) +
  theme_bw()

asv_trim %>%
  select(asv_code, seq_length) %>%
  left_join(asv_nseqs, by = "asv_code")

asv_trim %>%
  filter(seq_length < 210 | seq_length > 240) %>%
  select(asv_code, seq_length) %>%
  left_join(asv_nseqs, by = "asv_code") %>%
  summarize(nseqs = sum(nseqs),
            nasvs = n()) %>%
  mutate(pct_seqs = nseqs/sum(asv_melt$nseqs) * 100,
         pct_asvs = nasvs/length(unique(asv_melt$asv_code)) * 100)

# Remove ASVs with sequence lengths outside lower and upper limits
asv_rmpos <- asv_trim %>%
  filter(seq_length >= 210 & seq_length <= 240) %>%
  select(asv_code, trimmed) %>%
  arrange(asv_code)
length(unique(asv_rmpos$asv_code)) == length(unique(asv_rmpos$trimmed))  # Should evaluate to TRUE


#### Write new files for retained ASVs ####
# Write new fasta file for retained ASVs
# write.fasta(sequences = as.list(asv_rmpos$trimmed),
#             names = asv_rmpos$asv_code,
#             file.out = "output/fasta/lp2017-2019watercolumn18s_rmpos.fasta")

# Write new melted ASV table for retained ASVs
asv_melt_rmpos <- asv_melt %>%
  filter(asv_code %in% asv_rmpos$asv_code)
# asv_melt_rmpos %>%
#   write_tsv("output/melted/lp2017-2019watercolumn18s_melt_rmpos.tsv")
