# Format ENA read submission template

# Load libraries
library(tidyverse)
library(janitor)


#### Import and format data ####
# Import read information
reads <- read_tsv("data/ena/lp2017-2019watercolumn18s_read_files.txt", col_names = "file_name")

# Import read md5 information
md5 <- read_table("data/ena/lp2017-2019watercolumn18s.md5", col_names = c("md5", "file_name"))

# Import empty ENA template
ena <- read_tsv("data/ena/experiment_paired_fastq_spreadsheet.tsv", col_names = TRUE)


#### Format ENA template submission ####
# Format md5
md5 <- md5 %>%
  mutate(sample_alias = str_extract(file_name, "[:digit:][:digit:]\\-[:digit:][:digit:][:digit:]")) %>%
  mutate(sample_alias = str_c("LakePulse18S_water_", sample_alias)) %>%
  mutate(file_type = case_when(grepl("R1.fastq.gz$", file_name) ~ "forward_file_md5",
                               grepl("R2.fastq.gz$", file_name) ~ "reverse_file_md5")) %>%
  select(-file_name) %>%
  pivot_wider(names_from = file_type, values_from = md5)

# Relate read and md5 files to samples
samples <- reads %>%
  mutate(sample_alias = str_extract(file_name, "[:digit:][:digit:]\\-[:digit:][:digit:][:digit:]")) %>%
  mutate(sample_alias = str_c("LakePulse18S_water_", sample_alias)) %>%
  mutate(file_type = case_when(grepl("R1.fastq.gz$", file_name) ~ "forward_file_name",
                               grepl("R1.fastq.gz.md5$", file_name) ~ "forward_file_unencrypted_md5",
                               grepl("R2.fastq.gz$", file_name) ~ "reverse_file_name",
                               grepl("R2.fastq.gz.md5$", file_name) ~ "reverse_file_unencrypted_md5")) %>%
  pivot_wider(names_from = file_type, values_from = file_name) %>%
  left_join(md5, by = "sample_alias")

# Combine and complete fields
template <- ena %>%
  select(-sample_alias,
         -forward_file_name, -forward_file_md5, -forward_file_unencrypted_md5,
         -reverse_file_name, -reverse_file_md5, -reverse_file_unencrypted_md5) %>%
  bind_cols(samples) %>%
  select(project_accession, project_alias, sample_alias, experiment_alias,
         run_alias, library_name, library_source, library_selection,
         library_strategy, design_description, library_construction_protocol,
         instrument_model, file_type, library_layout, insert_size,
         forward_file_name, forward_file_md5, forward_file_unencrypted_md5,
         reverse_file_name, reverse_file_md5, reverse_file_unencrypted_md5) %>%
    mutate(instrument_model = "Illumina MiSeq",
           library_source = "METAGENOMIC",
           library_selection = "PCR",
           library_strategy = "AMPLICON",
           insert_size = 230)
  
# Write template data to file
# template %>%
#   write_tsv("output/ena/lp2017-2019watercolumn18s_ena_template_reads.tsv")
